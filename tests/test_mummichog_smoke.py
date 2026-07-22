"""Smoke test: a built model loads via `mummichog -n` and completes a run.

This is the authoritative acceptance check (MODEL_CONTRACT.md, Acceptance #1):
if mummichog's -n parses the JSON and the run exits 0 producing pathway + module
tables, the field names are right. Skipped when mummichog isn't installed.
"""

import json
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

pytest.importorskip("mass2chem", reason="scoped optional dep")
pytest.importorskip("metDataModel", reason="scoped optional dep")
pytest.importorskip("mummichog", reason="mummichog needed to validate the model")

from src import kegg_entities as ke  # noqa: E402
from src import prepare_mummichog_model as pm  # noqa: E402

# load_source's KEGG downloads live in kegg_entities now; patch them there.
FIX = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures")


def _build(tmp_path, monkeypatch):
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
    model_path, manifest_path = pm.prepare_mummichog_model(
        kegg_code="cre", out_dir=str(tmp_path), cache_dir=str(tmp_path),
        model_organism="Chlamydomonas reinhardtii",
        target_organism="Coelastrella sp.", date="20260711", validate=True)
    return model_path, manifest_path


def test_model_loads_and_runs_in_mummichog(tmp_path, monkeypatch):
    model_path, manifest_path = _build(tmp_path, monkeypatch)
    with open(manifest_path) as f:
        manifest = json.load(f)
    v = manifest["validation"]
    assert v["loads_via_-n"] is True
    assert v["smoke_run_exit_0"] is True
    assert v["mass_spotcheck_passed"] is True
    assert v["smoke_run"]["pathway_table_written"] is True
    assert v["smoke_run"]["module_table_written"] is True


def test_validate_model_helper_directly(tmp_path, monkeypatch):
    model_path, _ = _build(tmp_path, monkeypatch)
    result = pm.validate_model(model_path)
    assert result["smoke_run_exit_0"] is True
    assert result["pathway_table_written"] is True
