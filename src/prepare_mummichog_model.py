"""Build an organism-specific metabolic model for the mummichog pathway tool.

This is a NEW, self-contained output type. Unlike the gene-set modules
(pathway->gene / GO->gene / GMT) this model is COMPOUND-centric: compounds carry
a neutral formula + neutral monoisotopic mass, reactions carry real
substrate/product links, and pathways carry lists of reactions. It is pulled
directly from KEGG compound/reaction/pathway entities -- never derived from gene
sets -- and serialized to the metDataModel ``MetabolicModel`` shape that
``mummichog.main -n <model>.json`` consumes (mummichog 2.7.0).

Outputs (see MODEL_CONTRACT.md, the source of truth):

* ``<model_kegg_code>_<source>_<YYYYMMDD>.json``   -- the model
* ``<same-stem>.manifest.json``                    -- provenance, counts, sha256

Optionally (``emit_compound_sets=True``) the same ``load_source`` result is also
written as an ID-based compound-set GMT + readable table (see
:mod:`prepare_kegg_compound_sets`); those companions are recorded in the manifest
under ``companion_files``.

Input source is a config choice:

* ``kegg_org`` -- PRIMARY: all compounds/reactions/pathways for a KEGG organism
  code (e.g. ``cre`` for Chlamydomonas reinhardtii, used as a surrogate for
  Coelastrella which is not in KEGG). Implemented here.
* ``kaas``     -- SECONDARY: a non-model organism from a KAAS KO list
  (KO/EC -> reactions -> compounds -> pathways). A clean seam is left for it
  (see :func:`kegg_entities.load_source`); it is not implemented yet.

Heavy, model-specific dependencies (metDataModel for the structure, mass2chem
for the masses) are imported lazily so the gene-set modules stay dependency-free.
Pin them via ``requirements-mummichog.txt``. The source-loading step lives in
:mod:`kegg_entities` so the gene-set / compound-set modules can reuse it.
"""

import json
import os
import subprocess
import sys
from collections import defaultdict
from datetime import datetime

from .kegg_entities import SOURCE_TOKENS, is_metabolic_map, load_source
from .masses import neutral_mono_mass
from .utils import ensure_dir, git_head_sha, log_msg, sha256_file

PROTON = 1.00727646677  # matches mummichog.config.PROTON


# ---------------------------------------------------------------------------
# Assembly -- pure (no I/O). Uses metDataModel core_simple for the structure.
# ---------------------------------------------------------------------------

def assemble_model(reactions, compounds, pathway_names, *,
                   pathway_prefix, model_id, meta_data):
    """Assemble entities into a mummichog-ready model dict.

    Enforces the contract's hard requirements: every emitted compound has a
    neutral formula + mass (mummichog needs mass > 1); reactions reference only
    emitted compound ids; pathways reference only emitted reaction ids. Returns
    ``(model_dict, counts, spotcheck_pairs)`` where spotcheck_pairs is a list of
    (cid, computed_mass, kegg_exact_mass) for the mass spot-check.
    """
    # metDataModel core_simple is pandas-free and its field names match the
    # mummichog loader exactly; imported lazily as an optional dependency.
    try:
        from metDataModel.core_simple import (Compound, MetabolicModel, Pathway,
                                              Reaction)
    except ImportError as exc:  # pragma: no cover - exercised via install hint
        raise ImportError(
            "metDataModel is required for the mummichog model module. Install "
            "the scoped optional deps:  pip install -r requirements-mummichog.txt"
        ) from exc

    # 1. valid compounds: concrete neutral formula + monoisotopic mass > 1.
    valid, dropped_cpds = {}, 0
    for cid, info in compounds.items():
        mass = neutral_mono_mass(info.get("formula", ""))
        if mass is None or mass <= 1:
            dropped_cpds += 1
            continue
        c = Compound()
        c.id = cid
        c.name = info.get("name") or cid
        c.neutral_formula = info["formula"].strip()
        c.neutral_mono_mass = mass
        valid[cid] = c

    # 2. reactions pruned to valid compounds; keep only real substrate->product
    #    links (>=1 reactant AND >=1 product remaining).
    kept, used_cpds, rxn_refpaths, dropped_rxns = {}, set(), {}, 0
    for r in reactions:
        rr = [c for c in r["reactants"] if c in valid]
        pp = [c for c in r["products"] if c in valid]
        if not (rr and pp):
            dropped_rxns += 1
            continue
        obj = Reaction()
        obj.id = r["id"]
        obj.reactants = rr
        obj.products = pp
        obj.enzymes = list(r.get("enzymes", []))
        kept[r["id"]] = obj
        rxn_refpaths[r["id"]] = r.get("ref_pathways", [])
        used_cpds.update(rr)
        used_cpds.update(pp)

    # 3. emit only compounds actually used by a kept reaction (all connected).
    emit_cpds = [valid[c] for c in sorted(used_cpds)]

    # mass spot-check pairs (computed vs KEGG EXACT_MASS) for emitted compounds.
    spot = [(c.id, c.neutral_mono_mass, compounds[c.id]["exact_mass"])
            for c in emit_cpds if compounds[c.id].get("exact_mass") is not None]

    # 4. pathways: reference-pathway numbers -> organism pathway ids present in
    #    this organism's pathway list.
    pw_rxns = defaultdict(set)
    for rid, nums in rxn_refpaths.items():
        for num in nums:
            if not is_metabolic_map(num):  # drop global/overview maps
                continue
            pid = f"{pathway_prefix}{num}"
            if pid in pathway_names:
                pw_rxns[pid].add(rid)
    emit_pws = []
    for pid in sorted(pw_rxns):
        p = Pathway()
        p.id = pid
        p.name = pathway_names.get(pid, pid)
        p.list_of_reactions = sorted(pw_rxns[pid])
        emit_pws.append(p)

    # 5. metDataModel container -> contract JSON.
    mm = MetabolicModel()
    mm.id = model_id
    mm.meta_data = meta_data
    mm.list_of_compounds = emit_cpds
    mm.list_of_reactions = [kept[r] for r in sorted(kept)]
    mm.list_of_pathways = emit_pws

    counts = {
        "compounds": len(emit_cpds),
        "reactions": len(kept),
        "pathways": len(emit_pws),
        "compounds_dropped_no_mass": dropped_cpds,
        "reactions_dropped_unmapped": dropped_rxns,
    }
    return _emit(mm), counts, spot


def _emit(mm):
    """Serialize a metDataModel MetabolicModel to the exact mummichog JSON shape.

    Written explicitly (rather than via metDataModel's generic serialize) to
    keep the artifact clean, force ``identifiers`` to a dict (mummichog's
    fallback path calls ``.get`` on it) and surface ``meta_data.version`` (which
    the loader reads).
    """
    return {
        "id": mm.id,
        "version": mm.meta_data.get("version", ""),
        "meta_data": mm.meta_data,
        "list_of_compounds": [
            {"id": c.id, "name": c.name,
             "neutral_formula": c.neutral_formula,
             "neutral_mono_mass": c.neutral_mono_mass,
             "identifiers": {"kegg.compound": c.id}}
            for c in mm.list_of_compounds],
        "list_of_reactions": [
            {"id": r.id, "reactants": r.reactants, "products": r.products,
             "enzymes": r.enzymes}
            for r in mm.list_of_reactions],
        "list_of_pathways": [
            {"id": p.id, "name": p.name, "list_of_reactions": p.list_of_reactions}
            for p in mm.list_of_pathways],
    }


# ---------------------------------------------------------------------------
# Validation helpers
# ---------------------------------------------------------------------------

def evaluate_spotcheck(spot, tol_mda=1.0, min_n=5):
    """Compare our computed neutral masses to KEGG's own EXACT_MASS (independent).

    Passes iff >= *min_n* compounds were checked and all are within *tol_mda*
    milli-Da. KEGG rounds EXACT_MASS to 4 dp, well inside 1 mDa.
    """
    checks = sorted(
        ({"id": cid, "computed": computed, "kegg_exact_mass": kegg_exact,
          "abs_diff_mDa": round(abs(computed - kegg_exact) * 1000.0, 4)}
         for cid, computed, kegg_exact in spot),
        key=lambda d: d["id"])
    within = [c for c in checks if c["abs_diff_mDa"] <= tol_mda]
    return {
        "n_checked": len(checks),
        "n_within_1mDa": len(within),
        "max_abs_diff_mDa": max((c["abs_diff_mDa"] for c in checks), default=None),
        "tolerance_mDa": tol_mda,
        "passed": len(checks) >= min_n and len(within) == len(checks),
        "sample": checks[:10],
    }


def make_synthetic_feature_table(model, path, mode="pos", n_signif=15):
    """Write a tiny mummichog feature table (mz, rtime, p-value, t-score).

    Significant rows are the primary ion (M+H for pos, M-H for neg) of the
    model's own compounds, so at least one compound is guaranteed to be
    significant (mummichog's web report crashes on an all-empty significant
    list). Returns the number of significant features written.
    """
    cpds = [c for c in model["list_of_compounds"]
            if 50 <= c["neutral_mono_mass"] <= 2000][:n_signif]
    lines = ["mz\trtime\tp-value\tt-score"]
    rt = 60
    for i, c in enumerate(cpds):
        m = c["neutral_mono_mass"]
        mz = round(m - PROTON if mode == "neg" else m + PROTON, 5)
        lines.append(f"{mz}\t{rt}\t0.0001\t{4.0 + 0.1 * i:.2f}")
        rt += 4
    for j in range(12):  # background, clearly non-significant
        lines.append(f"{120.0 + j * 17.13:.5f}\t{50 + j * 6}\t{0.4 + 0.03 * j:.3f}\t0.5")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    return len(cpds)


def validate_model(model_path, mode="pos", permutations=20, cutoff=0.05):
    """Run ``mummichog -n <model>`` on a synthetic table derived from the model.

    Returns a dict for the manifest 'validation' block. Requires mummichog to be
    importable (raises ImportError otherwise). The acceptance criterion is that
    it loads via ``-n`` and completes (exit 0) producing pathway + module tables.
    """
    import importlib
    import shutil
    import tempfile

    version = getattr(importlib.import_module("mummichog.config"), "VERSION", None)
    with open(model_path) as f:
        model = json.load(f)
    workdir = tempfile.mkdtemp(prefix="mcg_validate_")
    try:
        n_sig = make_synthetic_feature_table(
            model, os.path.join(workdir, "features.txt"), mode=mode)
        cmd = [sys.executable, "-c",
               "from mummichog.command_line import main; main()",
               "-f", "features.txt", "-o", "validate",
               "-k", workdir, "-n", os.path.abspath(model_path),
               "-m", "neg" if mode == "neg" else "pos",
               "-u", "10", "-p", str(permutations), "-c", str(cutoff)]
        proc = subprocess.run(cmd, cwd=workdir, capture_output=True, text=True,
                              timeout=600)
        exit0 = proc.returncode == 0
        pathway_tbl = module_tbl = False
        for _root, _dirs, files in os.walk(workdir):
            for fn in files:
                if fn.endswith(".tsv") and "pathwayanalysis" in fn:
                    pathway_tbl = True
                if fn.endswith(".tsv") and "modularanalysis" in fn:
                    module_tbl = True
        return {
            "mummichog_version_tested": version,
            "loads_via_-n": exit0,
            "smoke_run_exit_0": exit0,
            "pathway_table_written": pathway_tbl,
            "module_table_written": module_tbl,
            "n_significant_features": n_sig,
            "stderr_tail": [] if exit0 else proc.stderr.strip().splitlines()[-6:],
        }
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


# ---------------------------------------------------------------------------
# Orchestration + output
# ---------------------------------------------------------------------------

def prepare_mummichog_model(kegg_code, out_dir, cache_dir, *,
                            source="kegg_org", model_organism=None,
                            model_kegg_code=None, target_organism=None,
                            model_is_surrogate=None, source_version=None,
                            date=None, refresh=False, validate=False,
                            kaas_file=None, emit_compound_sets=False):
    """Build + write the model JSON and its sidecar manifest.

    When *emit_compound_sets* is true, the same ``load_source`` result is also
    written as an ID-based compound-set GMT + readable table (same stem, same
    KEGG snapshot), recorded in the manifest under ``companion_files``.

    Returns ``(model_path, manifest_path)``.
    """
    date = date or datetime.utcnow().strftime("%Y%m%d")
    model_kegg_code = model_kegg_code or kegg_code
    if not model_kegg_code:
        raise ValueError(
            "model_kegg_code is required (a short label for the model file "
            "stem, e.g. 'cre'); for source='kaas' pass it explicitly")
    if model_is_surrogate is None:
        model_is_surrogate = bool(target_organism
                                  and target_organism != model_organism)
    stem = f"{model_kegg_code}_{SOURCE_TOKENS.get(source, source)}_{date}"

    src = load_source(source, kegg_code, cache_dir, refresh, kaas_file)
    resolved_version = source_version or src.get("source_version") or "unknown"

    # meta_data.version is REQUIRED by the mummichog loader (reads
    # meta_data['version']); the contract's meta_data block omits it, so we add
    # it here -- non-breaking, and MODEL_CONTRACT.md makes the loader the
    # authority on field names.
    meta_data = {
        "version": date,
        "model_organism": model_organism or "",
        "model_kegg_code": model_kegg_code,
        "target_organism": target_organism or "",
        "model_is_surrogate": model_is_surrogate,
        "source": src.get("source", "KEGG REST"),
        "source_version": resolved_version,
        "builder_version": git_head_sha() or "unknown",
    }

    model, counts, spot = assemble_model(
        src["reactions"], src["compounds"], src["pathway_names"],
        pathway_prefix=src["pathway_prefix"], model_id=stem, meta_data=meta_data)
    log_msg("model assembled: ", counts["compounds"], " compounds, ",
            counts["reactions"], " reactions, ", counts["pathways"], " pathways ",
            "(dropped ", counts["compounds_dropped_no_mass"], " compounds w/o mass, ",
            counts["reactions_dropped_unmapped"], " reactions unmapped)")

    ensure_dir(out_dir)
    model_path = os.path.join(out_dir, f"{stem}.json")
    with open(model_path, "w") as f:
        json.dump(model, f, indent=2)
    log_msg("wrote model: ", model_path)

    spotcheck = evaluate_spotcheck(spot)
    log_msg("mass spot-check: ", spotcheck["n_within_1mDa"], "/",
            spotcheck["n_checked"], " within 1 mDa (max ",
            spotcheck["max_abs_diff_mDa"], " mDa)")

    validation = {
        "mummichog_version_tested": "2.7.0",
        "loads_via_-n": None,
        "smoke_run_exit_0": None,
        "mass_spotcheck_passed": spotcheck["passed"],
    }
    if validate:
        try:
            v = validate_model(model_path)
            validation["mummichog_version_tested"] = v["mummichog_version_tested"]
            validation["loads_via_-n"] = v["loads_via_-n"]
            validation["smoke_run_exit_0"] = v["smoke_run_exit_0"]
            validation["smoke_run"] = v
            log_msg("mummichog -n validation: exit_0=", v["smoke_run_exit_0"],
                    " pathway_table=", v["pathway_table_written"])
        except ImportError as exc:
            log_msg("validation skipped (mummichog not installed): ", exc)
            validation["note"] = ("mummichog not installed; -n load/run not "
                                   "exercised (install mummichog==2.7.0 to enable)")

    manifest = {
        "model_file": os.path.basename(model_path),
        "sha256": sha256_file(model_path),
        "model_organism": meta_data["model_organism"],
        "model_kegg_code": model_kegg_code,
        "target_organism": meta_data["target_organism"],
        "model_is_surrogate": model_is_surrogate,
        "source": meta_data["source"],
        "source_version": resolved_version,
        "build_timestamp_utc": datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ"),
        "builder_git_sha": meta_data["builder_version"],
        "counts": {k: counts[k] for k in ("compounds", "reactions", "pathways")},
        "build_details": {
            "compounds_dropped_no_mass": counts["compounds_dropped_no_mass"],
            "reactions_dropped_unmapped": counts["reactions_dropped_unmapped"],
            "kegg_db_dates": src.get("kegg_db_dates", {}),
            "ko_coverage": src.get("ko_coverage", {}),
        },
        "mass_spotcheck": spotcheck,
        "validation": validation,
    }

    # Optional companion artifacts: ID-based compound-set GMT + readable table,
    # built from the SAME src (same KEGG snapshot). Recorded in the manifest so
    # the pipeline can pin/verify them alongside the model if it ever chooses to.
    if emit_compound_sets:
        from .prepare_kegg_compound_sets import (pathway_compound_sets,
                                                 write_compound_set_files)
        pw2cpd = pathway_compound_sets(
            src["reactions"], src["pathway_names"], src["compounds"],
            src["pathway_prefix"])
        written = write_compound_set_files(
            pw2cpd, src["pathway_names"], src["compounds"], out_dir, stem)
        manifest["companion_files"] = [
            {"file": os.path.basename(written["gmt_path"]),
             "role": "compound_pathway_gmt",
             "sha256": sha256_file(written["gmt_path"]),
             "counts": written["counts"]},
            {"file": os.path.basename(written["tab_path"]),
             "role": "pathway2compound_table",
             "sha256": sha256_file(written["tab_path"]),
             "counts": written["counts"]},
        ]

    manifest_path = os.path.join(out_dir, f"{stem}.manifest.json")
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)
    log_msg("wrote manifest: ", manifest_path)
    return model_path, manifest_path
