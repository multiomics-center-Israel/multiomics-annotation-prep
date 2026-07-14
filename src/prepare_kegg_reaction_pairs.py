"""Build a KEGG reaction-pair reference for the metabolite-network stage.

A NEW, self-contained output type. Unlike the mummichog model (organism-specific,
compound-centric) this reference is organism-INDEPENDENT: it is the whole KEGG
reaction database reduced to compound pairs, consumed by ``multiomic-core``'s
metabolite-network stage, which filters it to the compounds present in a study.

Pair definition -- ``equation_side_cartesian_product``: for each reaction, the
Cartesian product of the (de-duplicated) compound ids on the two sides of the
``<=>`` equation, dropping self-pairs. This reproduces the network the pipeline
drew before this reference existed. It is deliberately NOT KEGG "main reactant
pairs" and NOT RCLASS chemistry -- switching to those would change the biology
of the existing network and is out of scope here.

Outputs (all under *out_dir*, dated; nothing is published by this module):

* ``kegg_reaction_pairs_cross_side_<YYYYMMDD>.tsv.gz``       -- the reference
* ``kegg_reaction_pairs_cross_side_<YYYYMMDD>.manifest.json``-- provenance + sha256
* ``kegg_reaction_pairs_cross_side_<YYYYMMDD>.excluded.tsv`` -- per-reaction drops

Exclusion reason codes (reaction id + reason + diagnostic context):

* ``missing_equation``            -- no EQUATION field
* ``unsupported_equation_arrow``  -- EQUATION present, no recognized ``<=>`` arrow
* ``empty_left_compound_side``    -- no ``C#####`` compound left of the arrow
* ``empty_right_compound_side``   -- no ``C#####`` compound right of the arrow
* ``only_self_pairs``             -- every cross-pair was a self-pair (dropped)
* ``parse_error``                 -- an exception while turning the record into pairs
* ``fetch_failed``                -- record never parsed after retries (FATAL, see below)

Correctness guarantees:

* Reconciliation, not caching, defines retrieval success (see
  :func:`download_kegg.fetch_and_reconcile_reactions`).
* ``fetch_failed`` is FATAL by default: if any reaction id is unresolved the
  build refuses to finalize/publish an asset and writes only a failure report.
  ``allow_incomplete=True`` is a development-only escape that likewise produces
  NO publishable asset.
* The asset is assembled in a staging dir and atomically renamed into place only
  after the table, counts, manifest and sha256 are all complete, so a failed run
  never leaves a partial ``.tsv.gz``/``.manifest.json`` behind.
* Output rows are stably sorted, i.e. deterministically ORDERED for a fixed KEGG
  snapshot. Byte-identical compression is NOT claimed (``retrieved_at``, builder
  sha and gzip metadata vary between builds).
"""

import json
import os
import shutil
import tempfile
from collections import Counter
from datetime import datetime

from .download_kegg import (
    download_kegg_info,
    download_reaction_list,
    fetch_and_reconcile_reactions,
    kegg_snapshot_version,
    parse_kegg_info_dates,
    parse_reaction_list,
)
from .utils import ensure_dir, git_sha, log_msg, sha256_file, write_tab, write_tsv_gz

SCHEMA_VERSION = 1
PAIR_METHOD = "equation_side_cartesian_product"
STEM_PREFIX = "kegg_reaction_pairs_cross_side"
COLUMNS = ["reaction_id", "substrate_id", "product_id", "equation", "equation_arrow"]
EXCLUSION_COLUMNS = ["reaction_id", "reason", "context"]


def build_reaction_pairs(reactions):
    """Cartesian product of the two equation sides -> pair rows + exclusions.

    PURE (no I/O). For each parsed reaction record, emits one row per
    ``(reactant, product)`` pair -- both sides are already de-duplicated by
    :func:`download_kegg.parse_reaction_record`, so repeated ids within a side
    collapse -- dropping self-pairs (``substrate == product``). Compounds that
    appear on both sides keep their non-self cross-pairs. Returns
    ``(rows, excluded)``: rows are 5-tuples matching :data:`COLUMNS`; excluded
    are ``(reaction_id, reason, context)`` triples using the reason codes in the
    module docstring.
    """
    rows, excluded = [], []
    for r in reactions:
        rid = r.get("id", "")
        try:
            equation = r.get("equation", "") or ""
            arrow = r.get("equation_arrow")
            if not equation:
                excluded.append((rid, "missing_equation", ""))
                continue
            if not arrow:
                excluded.append((rid, "unsupported_equation_arrow", equation))
                continue
            reactants = r.get("reactants", [])
            products = r.get("products", [])
            if not reactants:
                excluded.append((rid, "empty_left_compound_side", equation))
                continue
            if not products:
                excluded.append((rid, "empty_right_compound_side", equation))
                continue
            pairs = [(rid, s, p, equation, arrow)
                     for s in reactants for p in products if s != p]
            if not pairs:
                excluded.append((rid, "only_self_pairs", equation))
                continue
            rows.extend(pairs)
        except Exception as exc:  # noqa: BLE001 - reaction kept as a reported drop
            excluded.append((rid, "parse_error", str(exc)))
    return rows, excluded


def _write_excluded(excluded, path):
    """Write the (reaction_id, reason, context) exclusion report."""
    write_tab(excluded, EXCLUSION_COLUMNS, path)


def _kegg_provenance(cache_dir, refresh):
    """Best-effort KEGG snapshot provenance (never a formal 'release' number).

    info/kegg carries only per-database dates; we record the snapshot string,
    the reaction database date, when we retrieved it, and the raw info lines.
    """
    info = {
        "source_snapshot": None,
        "reaction_database_date": None,
        "retrieved_at": datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ"),
        "kegg_info_raw": None,
    }
    try:
        info_path = download_kegg_info("kegg", cache_dir, refresh)
        with open(info_path) as f:
            info["kegg_info_raw"] = f.read().strip().splitlines()[:12]
        dates = parse_kegg_info_dates(info_path)
        info["source_snapshot"] = kegg_snapshot_version(dates)
        info["reaction_database_date"] = dates.get("reaction")
    except Exception as exc:  # noqa: BLE001 - provenance is best-effort
        log_msg("  (KEGG snapshot provenance unavailable: ", exc, ")")
    return info


def prepare_kegg_reaction_pairs(out_dir, cache_dir, *, date=None, refresh=False,
                                method=PAIR_METHOD, retries=3, rate_limit_s=0.34,
                                allow_incomplete=False):
    """Build + write the reaction-pair reference, its manifest and drop report.

    Returns ``(data_path, manifest_path, excluded_path)`` on a complete build.
    Raises RuntimeError (writing only the failure report) if any reaction id is
    ``fetch_failed`` and ``allow_incomplete`` is False. With ``allow_incomplete``
    it returns ``(None, None, report_path)`` and never writes a publishable asset.
    """
    date = date or datetime.utcnow().strftime("%Y%m%d")
    stem = f"{STEM_PREFIX}_{date}"
    ensure_dir(out_dir)

    # 1. every KEGG reaction id.
    reaction_ids = parse_reaction_list(download_reaction_list(cache_dir, refresh))
    log_msg("KEGG reactions to fetch: ", len(reaction_ids))

    # 2. fetch + reconcile -- retrieval success is decided by parsed ids, never
    #    by a bare cache hit (unresolved ids are retried cache-bypassed).
    parsed, missing = fetch_and_reconcile_reactions(
        reaction_ids, cache_dir, refresh=refresh, retries=retries,
        rate_limit_s=rate_limit_s)
    log_msg("reconciled: requested=", len(reaction_ids), " parsed=", len(parsed),
            " fetch_failed=", len(missing))

    # 3. build cross-side pairs (+ biological/parser exclusions).
    reactions = [parsed[rid] for rid in sorted(parsed)]
    rows, excluded = build_reaction_pairs(reactions)
    for rid in missing:
        excluded.append((rid, "fetch_failed", "no record parsed after retries"))

    # 4. deterministic ordering for a fixed snapshot.
    rows.sort(key=lambda t: (t[0], t[1], t[2]))
    excluded.sort(key=lambda t: (t[1], t[0]))

    # Reconciliation invariant: every requested id lands in exactly one bucket
    # (produced a pair, or was excluded for a documented reason / fetch_failed).
    produced_ids = {t[0] for t in rows}
    excluded_ids = {t[0] for t in excluded}
    unaccounted = set(reaction_ids) - produced_ids - excluded_ids
    if unaccounted:
        raise RuntimeError(
            f"{len(unaccounted)} reaction id(s) reconciled to neither a pair nor "
            f"an exclusion (e.g. {sorted(unaccounted)[:5]}); refusing to finalize.")

    by_reason = Counter(reason for _rid, reason, _ctx in excluded)

    # 5. fetch_failed is fatal by default; allow_incomplete stays NON-publishable.
    if missing:
        report_only = os.path.join(out_dir, f"{stem}.excluded.tsv")
        _write_excluded(excluded, report_only)
        msg = (f"{len(missing)} reaction(s) still fetch_failed after {retries} "
               f"retries; wrote failure report {report_only}.")
        if not allow_incomplete:
            raise RuntimeError(
                msg + " Refusing to finalize a reference asset -- re-run to resume "
                "(cached batches are reused; unresolved ids are re-fetched).")
        log_msg("WARNING (allow_incomplete, DEV ONLY): ", msg,
                " Producing NO publishable asset.")
        return None, None, report_only

    # 6. assemble in a staging dir; finalize by atomic rename only once the
    #    table, counts, manifest and sha256 are all complete.
    counts = {
        "n_rows": len(rows),
        "n_reactions": len(produced_ids),
        "n_compounds": len({c for _r, s, p, _e, _a in rows for c in (s, p)}),
    }
    prov = _kegg_provenance(cache_dir, refresh)
    builder = git_sha() or "unknown"

    staging = tempfile.mkdtemp(prefix=f".staging_{stem}_", dir=out_dir)
    try:
        data_stage = os.path.join(staging, f"{stem}.tsv.gz")
        excl_stage = os.path.join(staging, f"{stem}.excluded.tsv")
        manifest_stage = os.path.join(staging, f"{stem}.manifest.json")

        # Embedded header = the pipeline's canonical runtime metadata.
        header_lines = [
            f"schema_version: {SCHEMA_VERSION}",
            f"pair_definition_method: {method}",
            "source: KEGG REST",
            f"source_snapshot: {prov['source_snapshot']}",
            f"reaction_database_date: {prov['reaction_database_date']}",
            f"retrieved_at: {prov['retrieved_at']}",
            f"builder_git_sha: {builder}",
            f"n_rows: {counts['n_rows']}",
            f"n_reactions: {counts['n_reactions']}",
            f"n_compounds: {counts['n_compounds']}",
        ]
        write_tsv_gz(rows, COLUMNS, data_stage, header_lines=header_lines)
        _write_excluded(excluded, excl_stage)

        data_sha = sha256_file(data_stage)
        manifest = {
            "data_file": f"{stem}.tsv.gz",
            "sha256": data_sha,
            "schema_version": SCHEMA_VERSION,
            "pair_definition_method": method,
            "source": "KEGG REST",
            "source_snapshot": prov["source_snapshot"],
            "reaction_database_date": prov["reaction_database_date"],
            "retrieved_at": prov["retrieved_at"],
            "builder_git_sha": builder,
            "counts": counts,
            "exclusions": {
                "n_excluded_reactions": len(excluded),
                "by_reason": dict(sorted(by_reason.items())),
                "report_file": f"{stem}.excluded.tsv",
            },
            "kegg_info_raw": prov["kegg_info_raw"],
        }
        with open(manifest_stage, "w") as f:
            json.dump(manifest, f, indent=2)

        # Atomic publish (same filesystem -> os.replace is atomic; overwrites any
        # stale same-version artifact rather than leaving it in place).
        data_final = os.path.join(out_dir, f"{stem}.tsv.gz")
        excl_final = os.path.join(out_dir, f"{stem}.excluded.tsv")
        manifest_final = os.path.join(out_dir, f"{stem}.manifest.json")
        os.replace(data_stage, data_final)
        os.replace(excl_stage, excl_final)
        os.replace(manifest_stage, manifest_final)
    finally:
        shutil.rmtree(staging, ignore_errors=True)

    log_msg("wrote reference: ", data_final, "  (", counts["n_rows"], " rows, ",
            counts["n_reactions"], " reactions, ", counts["n_compounds"],
            " compounds; sha256 ", data_sha[:12], "...)")
    return data_final, manifest_final, excl_final
