"""End-to-end (offline) test of the model builder.

The KEGG download functions are monkeypatched to serve saved fixtures, so this
exercises the real orchestration + parsing + assembly + writers without any
network. Requires the scoped optional deps (metDataModel, mass2chem).
"""

import json
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

pytest.importorskip("mass2chem", reason="scoped optional dep")
pytest.importorskip("metDataModel", reason="scoped optional dep")

from src import prepare_mummichog_model as pm  # noqa: E402

FIX = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures")


@pytest.fixture
def built_model(tmp_path, monkeypatch):
    """Build the fixture model into tmp_path; return (model, manifest)."""
    rn_text = open(os.path.join(FIX, "kegg_reactions.txt")).read()
    cpd_text = open(os.path.join(FIX, "kegg_compounds.txt")).read()

    monkeypatch.setattr(pm, "download_kegg_org_ko_links",
                        lambda code, cache, refresh=False:
                        os.path.join(FIX, "link_ko_cre.txt"))
    monkeypatch.setattr(pm, "download_ko_reaction_links",
                        lambda cache, refresh=False:
                        os.path.join(FIX, "link_reaction_ko.txt"))
    monkeypatch.setattr(pm, "kegg_get_batched",
                        lambda prefix, ids, cache, refresh=False:
                        rn_text if prefix == "rn" else cpd_text)
    monkeypatch.setattr(pm, "download_kegg_org_pathways",
                        lambda code, cache, refresh=False:
                        os.path.join(FIX, "cre_pathways.txt" if code
                                     else "ref_pathways.txt"))
    monkeypatch.setattr(pm, "download_kegg_info",
                        lambda target, cache, refresh=False:
                        os.path.join(FIX, "info_kegg.txt"))

    model_path, manifest_path = pm.prepare_mummichog_model(
        kegg_code="cre", out_dir=str(tmp_path), cache_dir=str(tmp_path),
        source="kegg_org", model_organism="Chlamydomonas reinhardtii",
        target_organism="Coelastrella sp.", date="20260711")

    assert os.path.basename(model_path) == "cre_kegg_20260711.json"
    with open(model_path) as f:
        model = json.load(f)
    with open(manifest_path) as f:
        manifest = json.load(f)
    model["_path"] = model_path
    return model, manifest


def test_top_level_shape(built_model):
    model, _ = built_model
    for key in ("id", "version", "meta_data", "list_of_compounds",
                "list_of_reactions", "list_of_pathways"):
        assert key in model
    assert model["id"] == "cre_kegg_20260711"
    # meta_data.version is what the mummichog loader actually reads.
    assert model["meta_data"]["version"] == "20260711"


def test_meta_data_records_surrogate(built_model):
    model, _ = built_model
    md = model["meta_data"]
    assert md["model_organism"] == "Chlamydomonas reinhardtii"
    assert md["model_kegg_code"] == "cre"
    assert md["target_organism"] == "Coelastrella sp."
    assert md["model_is_surrogate"] is True
    assert md["source"] == "KEGG REST"
    assert md["source_version"] == "KEGG snapshot 2026-07-07"


def test_compounds_have_required_neutral_fields(built_model):
    model, _ = built_model
    assert model["list_of_compounds"]
    for c in model["list_of_compounds"]:
        assert c["neutral_formula"]
        assert isinstance(c["neutral_mono_mass"], float)
        assert c["neutral_mono_mass"] > 1
        assert c["identifiers"]["kegg.compound"] == c["id"]
    ids = {c["id"] for c in model["list_of_compounds"]}
    # only compounds used by a kept reaction are emitted
    assert ids == {"C00002", "C00008", "C00022", "C00031",
                   "C00074", "C00085", "C00668"}


def test_polymer_and_unmapped_compounds_dropped(built_model):
    model, _ = built_model
    ids = {c["id"] for c in model["list_of_compounds"]}
    assert "C00369" not in ids   # (C6H10O5)n polymer -> no computable mass
    assert "C05345" not in ids   # never present as a compound


def test_reactions_pruned_to_real_links(built_model):
    model, _ = built_model
    rids = {r["id"] for r in model["list_of_reactions"]}
    assert rids == {"R00299", "R00771", "R00200"}
    # R01786 dropped: its only product C05345 is unmapped
    # R02110 dropped: its only reactant C00369 is a dropped polymer
    assert "R01786" not in rids
    assert "R02110" not in rids


def test_reactions_reference_existing_compounds(built_model):
    model, _ = built_model
    cpd_ids = {c["id"] for c in model["list_of_compounds"]}
    for r in model["list_of_reactions"]:
        assert r["reactants"] and r["products"]
        for c in r["reactants"] + r["products"]:
            assert c in cpd_ids


def test_pathways_reference_existing_reactions(built_model):
    model, _ = built_model
    rxn_ids = {r["id"] for r in model["list_of_reactions"]}
    pw = {p["id"]: p for p in model["list_of_pathways"]}
    assert set(pw) == {"cre00010"}
    assert pw["cre00010"]["name"] == "Glycolysis / Gluconeogenesis"
    assert set(pw["cre00010"]["list_of_reactions"]) == {"R00299", "R00771", "R00200"}
    for p in model["list_of_pathways"]:
        for rid in p["list_of_reactions"]:
            assert rid in rxn_ids


# KEGG "Global and overview maps" (01100-01299) that must never appear.
OVERVIEW_MAPS = ["01100", "01110", "01200", "01210", "01212", "01230",
                 "01232", "01240", "01250"]


def test_global_overview_maps_excluded(built_model):
    model, _ = built_model
    pids = {p["id"] for p in model["list_of_pathways"]}
    # the fixtures deliberately include cre01100 / cre01200 in the pathway list
    # and on reaction records; none of the overview maps may survive
    for num in OVERVIEW_MAPS:
        assert f"cre{num}" not in pids
    assert all(int(p["id"][3:]) < 1000 for p in model["list_of_pathways"])


def test_manifest_counts_and_checksum(built_model):
    import hashlib
    model, manifest = built_model
    assert manifest["counts"]["compounds"] == len(model["list_of_compounds"])
    assert manifest["counts"]["reactions"] == len(model["list_of_reactions"])
    assert manifest["counts"]["pathways"] == len(model["list_of_pathways"])
    assert manifest["model_is_surrogate"] is True
    # sha256 in the manifest matches the model file on disk
    with open(model["_path"], "rb") as f:
        assert manifest["sha256"] == hashlib.sha256(f.read()).hexdigest()


def test_mass_spotcheck_passed(built_model):
    _, manifest = built_model
    sc = manifest["mass_spotcheck"]
    assert sc["n_checked"] >= 5
    assert sc["passed"] is True
    assert sc["max_abs_diff_mDa"] <= 1.0


def test_manifest_records_ko_coverage(built_model):
    # Organism reactions are resolved gene->KO->reaction; coverage is reported.
    _, manifest = built_model
    cov = manifest["build_details"]["ko_coverage"]
    assert cov["n_genes"] == 6
    assert cov["n_kos"] == 5                 # K00844 de-duplicated across 2 genes
    assert cov["n_kos_with_reaction"] == 4   # K99999 maps to no reaction
    assert cov["n_reactions_from_kos"] == 5  # R00299,R00771,R01786,R00200,R02110


def test_manifest_records_kegg_snapshot_dates(built_model):
    _, manifest = built_model
    dates = manifest["build_details"]["kegg_db_dates"]
    assert dates == {"pathway": "2026-07-03", "reaction": "2026-07-07",
                     "compound": "2026-07-06"}


def test_kaas_source_builds_from_ko_file(tmp_path, monkeypatch):
    """The KAAS path is the same KO-list pipeline, KOs read from a file and
    pathways taken from the KEGG-wide reference maps (organism not in KEGG)."""
    rn_text = open(os.path.join(FIX, "kegg_reactions.txt")).read()
    cpd_text = open(os.path.join(FIX, "kegg_compounds.txt")).read()
    monkeypatch.setattr(pm, "download_ko_reaction_links",
                        lambda cache, refresh=False:
                        os.path.join(FIX, "link_reaction_ko.txt"))
    monkeypatch.setattr(pm, "kegg_get_batched",
                        lambda prefix, ids, cache, refresh=False:
                        rn_text if prefix == "rn" else cpd_text)
    monkeypatch.setattr(pm, "download_kegg_org_pathways",
                        lambda code, cache, refresh=False:
                        os.path.join(FIX, "cre_pathways.txt" if code
                                     else "ref_pathways.txt"))
    monkeypatch.setattr(pm, "download_kegg_info",
                        lambda target, cache, refresh=False:
                        os.path.join(FIX, "info_kegg.txt"))

    model_path, manifest_path = pm.prepare_mummichog_model(
        kegg_code=None, out_dir=str(tmp_path), cache_dir=str(tmp_path),
        source="kaas", kaas_file=os.path.join(FIX, "kaas_query.ko.txt"),
        model_organism="Coelastrella sp.", model_kegg_code="coel",
        date="20260711")

    assert os.path.basename(model_path) == "coel_kaas_20260711.json"
    with open(model_path) as f:
        model = json.load(f)
    with open(manifest_path) as f:
        manifest = json.load(f)
    # reference (map#####) pathways, since a non-model organism has none of its
    # own; the overview map (map01100/map01200) is filtered out here too
    assert {p["id"] for p in model["list_of_pathways"]} == {"map00010"}
    assert {r["id"] for r in model["list_of_reactions"]} == {"R00299", "R00771", "R00200"}
    assert model["meta_data"]["source"] == "KEGG REST (KAAS KO list)"
    assert model["meta_data"]["source_version"] == "KEGG snapshot 2026-07-07"
    cov = manifest["build_details"]["ko_coverage"]
    assert cov["n_kos"] == 4                 # K99999 has no reaction
    assert cov["n_kos_with_reaction"] == 3
    assert cov["n_reactions_from_kos"] == 4  # R00299,R00771,R01786,R00200
