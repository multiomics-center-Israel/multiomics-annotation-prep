"""KEGG compound-set GMT + readable table for ID-based metabolomics enrichment.

`multiomic-core` runs metabolomics enrichment two ways for a non-model organism:
m/z-based (mummichog ``-n <model>.json``, built by :mod:`prepare_mummichog_model`)
and ID-based (ORA / GSEA / QEA over a GMT of pathway compound sets). This module
produces the ID-based inputs from the SAME organism source used for the model
(:func:`kegg_entities.load_source`), so both enrichment paths describe the same
pathways and the same biology.

Two files are written (stem matches the model's:
``<model_kegg_code>_<source_token>_<YYYYMMDD>``):

* ``<stem>.compound_pathway.gmt``   -- ``pathway<tab>name<tab>C-ID...`` (GMT),
  consumed by multiomic-core's ``read_gmt`` / ``load_gene_sets``.
* ``<stem>.pathway2compound.tab``   -- readable table
  (pathway_id, pathway_name, compound_id, compound_name).

The compound set is a SUPERSET of the mummichog model's compounds: it is built
from the reactions' substrate/product links with NO mass filter (ID-based
enrichment needs no mass), whereas the model drops compounds without a
computable neutral mass. Same pathways, same biology, broader compound coverage.

Dependency-light: stdlib + the ``requests``-backed downloads via
:mod:`kegg_entities` / :mod:`download_kegg`. No ``metDataModel`` / ``mass2chem``.
"""

import json
import os
from collections import defaultdict
from datetime import datetime

from .kegg_entities import SOURCE_TOKENS, is_metabolic_map, load_source
from .utils import (ensure_dir, git_head_sha, log_msg, sha256_file, write_gmt,
                    write_pathway2compound)


# ---------------------------------------------------------------------------
# Pure core (no I/O, no network) -- easy to unit-test with synthetic records.
# ---------------------------------------------------------------------------

def pathway_compound_sets(reactions, pathway_names, compounds, pathway_prefix):
    """Build ``{pathway_id: set(compound_id)}`` from reaction substrate/products.

    For each reaction, for each of its reference-pathway numbers that is a real
    metabolic map (:func:`kegg_entities.is_metabolic_map`) AND resolves to a
    pathway present for this organism (``f"{pathway_prefix}{num}"`` in
    ``pathway_names``), the reaction's reactants + products are added to that
    pathway's compound set. Mirrors ``assemble_model``'s pathway construction but
    collects COMPOUNDS and applies NO mass filter, so the sets are a superset of
    the mummichog model's compounds (correct for ID-based ORA/GSEA/QEA).

    ``compounds`` is unused here (membership is driven by the reactions); it is
    accepted for signature symmetry with the writers and future variants.
    """
    pw2cpd = defaultdict(set)
    for r in reactions:
        members = r.get("reactants", []) + r.get("products", [])
        if not members:
            continue
        for num in r.get("ref_pathways", []):
            if not is_metabolic_map(num):
                continue
            pid = f"{pathway_prefix}{num}"
            if pid in pathway_names:
                pw2cpd[pid].update(members)
    return dict(pw2cpd)


def _set_counts(pw2cpd):
    """Summary counts for a pathway->compounds mapping (non-empty sets only)."""
    return {
        "pathways": sum(1 for v in pw2cpd.values() if v),
        "compounds": len({c for v in pw2cpd.values() for c in v}),
        "pairs": sum(len(v) for v in pw2cpd.values()),
    }


# ---------------------------------------------------------------------------
# File writing (shared by the standalone entry and prepare_mummichog_model's
# emit_compound_sets path -- one code path, one snapshot).
# ---------------------------------------------------------------------------

def write_compound_set_files(pw2cpd, pathway_names, compounds, out_dir, stem):
    """Write ``<stem>.compound_pathway.gmt`` + ``<stem>.pathway2compound.tab``.

    Returns ``{"gmt_path", "tab_path", "counts"}``.
    """
    ensure_dir(out_dir)
    gmt_path = os.path.join(out_dir, f"{stem}.compound_pathway.gmt")
    tab_path = os.path.join(out_dir, f"{stem}.pathway2compound.tab")

    # write_gmt de-dups + sorts members and skips empty sets.
    write_gmt(pw2cpd, pathway_names, gmt_path)

    rows = []
    for pid in sorted(pw2cpd):
        pname = pathway_names.get(pid, "")
        for cid in sorted(pw2cpd[pid]):
            cname = (compounds.get(cid, {}) or {}).get("name", "") if compounds else ""
            rows.append((pid, pname, cid, cname or ""))
    write_pathway2compound(rows, tab_path)

    return {"gmt_path": gmt_path, "tab_path": tab_path,
            "counts": _set_counts(pw2cpd)}


# ---------------------------------------------------------------------------
# Standalone orchestration (fetch -> build -> write -> sidecar manifest).
# ---------------------------------------------------------------------------

def prepare_kegg_compound_sets(kegg_code, out_dir, cache_dir, *,
                               source="kegg_org", model_kegg_code=None,
                               model_organism=None, target_organism=None,
                               source_version=None, date=None, refresh=False,
                               kaas_file=None):
    """Build + write the compound-set GMT, readable table, and sidecar manifest.

    Returns ``(gmt_path, tab_path, manifest_path)``.
    """
    date = date or datetime.utcnow().strftime("%Y%m%d")
    model_kegg_code = model_kegg_code or kegg_code
    if not model_kegg_code:
        raise ValueError(
            "model_kegg_code is required (a short label for the file stem, "
            "e.g. 'cre'); for source='kaas' pass it explicitly")
    stem = f"{model_kegg_code}_{SOURCE_TOKENS.get(source, source)}_{date}"

    src = load_source(source, kegg_code, cache_dir, refresh, kaas_file)
    resolved_version = source_version or src.get("source_version") or "unknown"

    pw2cpd = pathway_compound_sets(
        src["reactions"], src["pathway_names"], src["compounds"],
        src["pathway_prefix"])
    written = write_compound_set_files(
        pw2cpd, src["pathway_names"], src["compounds"], out_dir, stem)
    counts = written["counts"]
    log_msg("compound sets: ", counts["pathways"], " pathways, ",
            counts["compounds"], " compounds, ", counts["pairs"], " pairs")

    gmt_name = os.path.basename(written["gmt_path"])
    tab_name = os.path.basename(written["tab_path"])
    manifest = {
        "compound_pathway_gmt": gmt_name,
        "pathway2compound_table": tab_name,
        "sha256": {
            gmt_name: sha256_file(written["gmt_path"]),
            tab_name: sha256_file(written["tab_path"]),
        },
        "model_organism": model_organism or "",
        "model_kegg_code": model_kegg_code,
        "target_organism": target_organism or "",
        "model_is_surrogate": bool(target_organism
                                   and target_organism != model_organism),
        "source": src.get("source", "KEGG REST"),
        "source_version": resolved_version,
        "build_timestamp_utc": datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ"),
        "builder_git_sha": git_head_sha() or "unknown",
        "counts": counts,
        "build_details": {
            "kegg_db_dates": src.get("kegg_db_dates", {}),
            "ko_coverage": src.get("ko_coverage", {}),
        },
    }
    manifest_path = os.path.join(out_dir, f"{stem}.compound_sets.manifest.json")
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)
    log_msg("wrote manifest: ", manifest_path)
    return written["gmt_path"], written["tab_path"], manifest_path
