"""Tests for the KEGG compound-set GMT + table (prepare_kegg_compound_sets).

The pure-function tests use tiny synthetic in-memory records -- no network, and
NO mass/metDataModel deps (that's the whole point: ID-based enrichment needs no
mass). The end-to-end test reuses the same saved fixtures as the model tests,
monkeypatching the KEGG downloads on `kegg_entities` (where load_source resolves
them), and shows the compound set is a SUPERSET of the mummichog model's
compounds for the same pathway.
"""

import hashlib
import json
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src import kegg_entities as ke  # noqa: E402
from src import prepare_kegg_compound_sets as cs  # noqa: E402

FIX = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures")


# --- synthetic records for the pure functions --------------------------------

# 00010 -> present in names (org00010); 00500 -> metabolic but absent from names;
# 01200 -> overview map (>=1000), never a pathway.
SYNTH_REACTIONS = [
    {"id": "R1", "reactants": ["C1", "C2"], "products": ["C3"],
     "enzymes": [], "ref_pathways": ["00010", "01100"]},
    {"id": "R2", "reactants": ["C3"], "products": ["C4"],
     "enzymes": [], "ref_pathways": ["00010"]},
    {"id": "R3", "reactants": ["C5"], "products": ["C6"],
     "enzymes": [], "ref_pathways": ["00500"]},   # 00500 not in names -> dropped
    {"id": "R4", "reactants": ["Cx"], "products": ["Cy"],
     "enzymes": [], "ref_pathways": ["01200"]},    # overview map -> dropped
]
SYNTH_PATHWAY_NAMES = {"org00010": "Glycolysis", "org00020": "TCA cycle"}
SYNTH_COMPOUNDS = {
    "C1": {"id": "C1", "name": "Alpha", "formula": "C6H12O6", "exact_mass": 180.06},
    "C2": {"id": "C2", "name": "Beta", "formula": "", "exact_mass": None},  # massless
    # C3/C4 intentionally absent from the compounds dict -> name falls back to ""
}
PREFIX = "org"


def test_groups_compounds_by_pathway():
    pw2cpd = cs.pathway_compound_sets(
        SYNTH_REACTIONS, SYNTH_PATHWAY_NAMES, SYNTH_COMPOUNDS, PREFIX)
    assert set(pw2cpd) == {"org00010"}
    assert pw2cpd["org00010"] == {"C1", "C2", "C3", "C4"}


def test_overview_maps_and_absent_pathways_excluded():
    pw2cpd = cs.pathway_compound_sets(
        SYNTH_REACTIONS, SYNTH_PATHWAY_NAMES, SYNTH_COMPOUNDS, PREFIX)
    assert "org01200" not in pw2cpd   # R4: global/overview map (>=1000)
    assert "org00500" not in pw2cpd   # R3: metabolic but absent from names
    assert "org00020" not in pw2cpd   # no reaction maps here -> no empty set


def test_massless_compound_kept_superset_behavior():
    # C2 has no formula, so the mummichog model would drop it; the ID-based set
    # keeps it (no mass filter). This is the defining superset property.
    pw2cpd = cs.pathway_compound_sets(
        SYNTH_REACTIONS, SYNTH_PATHWAY_NAMES, SYNTH_COMPOUNDS, PREFIX)
    assert "C2" in pw2cpd["org00010"]


def test_write_files_format(tmp_path):
    pw2cpd = cs.pathway_compound_sets(
        SYNTH_REACTIONS, SYNTH_PATHWAY_NAMES, SYNTH_COMPOUNDS, PREFIX)
    written = cs.write_compound_set_files(
        pw2cpd, SYNTH_PATHWAY_NAMES, SYNTH_COMPOUNDS, str(tmp_path), "org_kegg_20260101")

    # GMT: TERM<tab>name<tab>member... , members sorted, empty sets skipped.
    with open(written["gmt_path"]) as f:
        gmt_lines = [ln.rstrip("\n") for ln in f if ln.strip()]
    assert gmt_lines == ["org00010\tGlycolysis\tC1\tC2\tC3\tC4"]

    # readable table: header + sorted rows, name falls back to "" when unknown.
    with open(written["tab_path"]) as f:
        tab_lines = [ln.rstrip("\n") for ln in f]
    assert tab_lines[0] == "pathway_id\tpathway_name\tcompound_id\tcompound_name"
    assert tab_lines[1] == "org00010\tGlycolysis\tC1\tAlpha"
    assert tab_lines[2] == "org00010\tGlycolysis\tC2\tBeta"
    assert tab_lines[3] == "org00010\tGlycolysis\tC3\t"   # C3 absent -> blank name
    assert tab_lines[4] == "org00010\tGlycolysis\tC4\t"

    assert written["counts"] == {"pathways": 1, "compounds": 4, "pairs": 4}


# --- end-to-end from the saved fixtures (offline, no mass deps) ---------------

def _patch_downloads(monkeypatch):
    rn_text = open(os.path.join(FIX, "kegg_reactions.txt")).read()
    cpd_text = open(os.path.join(FIX, "kegg_compounds.txt")).read()
    monkeypatch.setattr(ke, "download_kegg_org_ko_links",
                        lambda code, cache, refresh=False:
                        os.path.join(FIX, "link_ko_cre.txt"))
    monkeypatch.setattr(ke, "download_ko_reaction_links",
                        lambda cache, refresh=False:
                        os.path.join(FIX, "link_reaction_ko.txt"))
    monkeypatch.setattr(ke, "kegg_get_batched",
                        lambda prefix, ids, cache, refresh=False:
                        rn_text if prefix == "rn" else cpd_text)
    monkeypatch.setattr(ke, "download_kegg_org_pathways",
                        lambda code, cache, refresh=False:
                        os.path.join(FIX, "cre_pathways.txt" if code
                                     else "ref_pathways.txt"))
    monkeypatch.setattr(ke, "download_kegg_info",
                        lambda target, cache, refresh=False:
                        os.path.join(FIX, "info_kegg.txt"))


# the mummichog model emits these compounds for cre00010 (see test_model_assembly)
MODEL_CRE00010_COMPOUNDS = {"C00002", "C00008", "C00022", "C00031",
                            "C00074", "C00085", "C00668"}


def test_end_to_end_from_fixtures(tmp_path, monkeypatch):
    _patch_downloads(monkeypatch)
    gmt_path, tab_path, manifest_path = cs.prepare_kegg_compound_sets(
        kegg_code="cre", out_dir=str(tmp_path), cache_dir=str(tmp_path),
        source="kegg_org", model_organism="Chlamydomonas reinhardtii",
        target_organism="Coelastrella sp.", date="20260711")

    # stem matches the model's, with a content-descriptor suffix
    assert os.path.basename(gmt_path) == "cre_kegg_20260711.compound_pathway.gmt"
    assert os.path.basename(tab_path) == "cre_kegg_20260711.pathway2compound.tab"

    sets = {}
    with open(gmt_path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            sets[parts[0]] = (parts[1], set(parts[2:]))

    # only the organism metabolic map that has reactions survives; the overview
    # maps (cre01100/cre01200) and pathways absent from the org list are dropped
    assert set(sets) == {"cre00010"}
    name, members = sets["cre00010"]
    assert name == "Glycolysis / Gluconeogenesis"
    assert members == MODEL_CRE00010_COMPOUNDS | {"C05345"}

    # SUPERSET of the model: the extra compound is C05345, kept because the
    # ID-based set applies no mass filter (R01786 -> C05345 is not pruned)
    assert MODEL_CRE00010_COMPOUNDS < members
    assert members - MODEL_CRE00010_COMPOUNDS == {"C05345"}

    # readable table: header + one row per (pathway, compound)
    with open(tab_path) as f:
        tab_lines = [ln.rstrip("\n") for ln in f if ln.strip()]
    assert tab_lines[0] == "pathway_id\tpathway_name\tcompound_id\tcompound_name"
    assert len(tab_lines) - 1 == len(members)   # 8 data rows

    # sidecar manifest: counts, surrogate flag, and a sha256 per file
    with open(manifest_path) as f:
        manifest = json.load(f)
    assert manifest["counts"] == {"pathways": 1, "compounds": 8, "pairs": 8}
    assert manifest["model_is_surrogate"] is True
    assert set(manifest["sha256"]) == {os.path.basename(gmt_path),
                                       os.path.basename(tab_path)}


def test_emit_compound_sets_from_model_builder(tmp_path, monkeypatch):
    """Hybrid path: one prepare_mummichog_model run writes the model AND the
    compound-set companions from the same load_source, recorded in the model
    manifest's companion_files. Needs the model's scoped optional deps."""
    pytest.importorskip("mass2chem", reason="scoped optional dep")
    pytest.importorskip("metDataModel", reason="scoped optional dep")
    from src import prepare_mummichog_model as pm

    _patch_downloads(monkeypatch)
    model_path, manifest_path = pm.prepare_mummichog_model(
        kegg_code="cre", out_dir=str(tmp_path), cache_dir=str(tmp_path),
        source="kegg_org", model_organism="Chlamydomonas reinhardtii",
        target_organism="Coelastrella sp.", date="20260711",
        emit_compound_sets=True)

    gmt = tmp_path / "cre_kegg_20260711.compound_pathway.gmt"
    tab = tmp_path / "cre_kegg_20260711.pathway2compound.tab"
    assert gmt.exists() and tab.exists()

    with open(manifest_path) as f:
        manifest = json.load(f)
    comp = {c["role"]: c for c in manifest["companion_files"]}
    assert set(comp) == {"compound_pathway_gmt", "pathway2compound_table"}
    # each companion's recorded sha256 matches the file on disk
    for role, entry in comp.items():
        p = tmp_path / entry["file"]
        assert entry["sha256"] == hashlib.sha256(p.read_bytes()).hexdigest()
    # same compound universe as the standalone build (8 compounds for cre00010)
    assert comp["compound_pathway_gmt"]["counts"]["compounds"] == 8
