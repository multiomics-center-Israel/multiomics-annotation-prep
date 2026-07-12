"""Tests for the MODELS.md registry append logic (src/model_registry.py).

Pure text/JSON: mock manifests, no network and no optional deps.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.model_registry import (  # noqa: E402
    DEFAULT_REPO,
    append_model_row,
    build_row,
    existing_tags,
    tag_from_manifest,
)

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

CRE_MANIFEST = {
    "model_file": "cre_kegg_20260711.json",
    "sha256": "c403c96fbec8df9ae34b828fec01270c8ea3940acc36e4e5ff770868dc8b912b",
    "model_organism": "Chlamydomonas reinhardtii",
    "model_kegg_code": "cre",
    "target_organism": "Coelastrella sp.",
    "model_is_surrogate": True,
    "source": "KEGG REST",
    "source_version": "KEGG snapshot 2026-07-10",
    "build_timestamp_utc": "2026-07-11T09:00:00Z",
    "build_details": {"kegg_db_dates": {"pathway": "2026-07-03",
                                        "reaction": "2026-07-07",
                                        "compound": "2026-07-06"}},
    "counts": {"compounds": 1110, "reactions": 900, "pathways": 80},
}


def test_build_row_fields():
    row = build_row(CRE_MANIFEST)
    assert row == [
        "Coelastrella sp.",
        "Chlamydomonas reinhardtii (cre)",
        "yes",
        "cre_kegg_20260711",
        "https://github.com/multiomics-center-Israel/multiomics-annotation-prep/"
        "releases/download/cre_kegg_20260711/cre_kegg_20260711.json",
        "c403c96fbec8df9ae34b828fec01270c8ea3940acc36e4e5ff770868dc8b912b",
        "2026-07-10",   # snapshot, parsed from source_version
        "2026-07-11",   # build date, from the tag's YYYYMMDD
    ]


def test_append_creates_file_and_row(tmp_path):
    models = str(tmp_path / "MODELS.md")
    assert append_model_row(models, CRE_MANIFEST) is True
    text = open(models).read()
    assert "| Target organism |" in text          # header written
    assert "cre_kegg_20260711" in text
    assert "Coelastrella sp." in text
    assert existing_tags(text) == {"cre_kegg_20260711"}


def test_append_is_idempotent_on_tag(tmp_path):
    models = str(tmp_path / "MODELS.md")
    assert append_model_row(models, CRE_MANIFEST) is True
    assert append_model_row(models, CRE_MANIFEST) is False   # same tag -> skip
    text = open(models).read()
    rows = [ln for ln in text.splitlines()
            if ln.lstrip().startswith("|") and "cre_kegg_20260711" in ln]
    assert len(rows) == 1   # exactly one row for the tag


def test_append_second_distinct_model(tmp_path):
    models = str(tmp_path / "MODELS.md")
    append_model_row(models, CRE_MANIFEST)
    other = dict(CRE_MANIFEST, model_file="mmu_kegg_20260801.json",
                 model_kegg_code="mmu", model_organism="Mus musculus",
                 target_organism="", model_is_surrogate=False,
                 build_timestamp_utc="2026-08-01T00:00:00Z")
    assert append_model_row(models, other) is True
    tags = existing_tags(open(models).read())
    assert tags == {"cre_kegg_20260711", "mmu_kegg_20260801"}


def test_non_surrogate_target_falls_back_to_model_organism():
    m = dict(CRE_MANIFEST, target_organism="", model_is_surrogate=False,
             model_organism="Mus musculus", model_kegg_code="mmu")
    row = build_row(m)
    assert row[0] == "Mus musculus"   # target falls back to the model organism
    assert row[2] == "no"


def test_snapshot_falls_back_to_db_dates_when_no_source_version():
    m = dict(CRE_MANIFEST, source_version="")
    # newest of pathway/reaction/compound
    assert build_row(m)[6] == "2026-07-07"


def test_seed_row_matches_code_output():
    """The committed MODELS.md seed row must equal what build_row() produces,
    so publishing cre again is a genuine no-op (guards format drift)."""
    text = open(os.path.join(REPO_ROOT, "MODELS.md")).read()
    assert "cre_kegg_20260711" in existing_tags(text)
    expected = build_row(CRE_MANIFEST)
    seed = [ln for ln in text.splitlines() if "cre_kegg_20260711" in ln
            and ln.lstrip().startswith("|")]
    assert len(seed) == 1
    cells = [c.strip() for c in seed[0].strip().strip("|").split("|")]
    assert cells == [str(c) for c in expected]


def test_publishing_seeded_tag_is_noop(tmp_path):
    """Re-publishing the already-listed cre tag adds nothing."""
    models = str(tmp_path / "MODELS.md")
    # start from a copy of the real registry
    with open(os.path.join(REPO_ROOT, "MODELS.md")) as f:
        open(models, "w").write(f.read())
    assert append_model_row(models, CRE_MANIFEST, DEFAULT_REPO) is False
