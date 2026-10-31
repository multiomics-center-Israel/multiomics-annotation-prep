"""Tests for src/annotation_manifest.py (manifest.json in annotation_dir releases).

Offline: the cache dir is pre-seeded with go-basic.obo (the build's own file);
the live version lookups go through a fake cached_download that records each
call, and the "unavailable" case replaces it with one that raises.
"""

import json
import os
import subprocess
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from src import annotation_manifest as am  # noqa: E402
from src.download_kegg import parse_kegg_info_dates  # noqa: E402
from src.utils import sha256_file  # noqa: E402

FIXTURES = os.path.join(ROOT, "tests", "fixtures")

_OBO = "format-version: 1.2\ndata-version: releases/2026-01-01\n\n[Term]\nid: GO:0000001\n"

_REGISTRY = """<MartRegistry>
<MartURLLocation database="ensembl_mart_113" default="1" displayName="Ensembl Genes 113"
  host="www.ensembl.org" name="ENSEMBL_MART_ENSEMBL" path="/biomart/martservice"
  port="443" serverVirtualSchema="default" visible="1" />
<MartURLLocation database="plants_mart_60" displayName="Ensembl Plants Genes 60"
  host="plants.ensembl.org" name="plants_mart" serverVirtualSchema="plants_mart" visible="1" />
</MartRegistry>
"""

# A GAF-style build plus KEGG: one file of every kind the builders write.
_BUILT = {
    "GO2gene_BP.tab": "GO\tGene\nGO:0000001\tg1\n",
    "GO2name_BP.tab": "GO\tTerm\nGO:0000001\troot process\n",
    "GO_BP.gmt": "GO:0000001\troot process\tg1\n",
    "Annotation.tab": "Gene\tgene_name\tprotein_name\ng1\tABC1\tprotein\n",
    "KEGG_pathway2gene.tab": "v1\tindex\nmap00010\tg1\n",
    "KEGG_pathway2name.tab": "pathway\tinfo\nmap00010\tGlycolysis\n",
    "GO_provenance.json": "{}\n",
}


def _write(path, text):
    with open(path, "w") as f:
        f.write(text)


@pytest.fixture
def build(tmp_path):
    out = tmp_path / "annotation_out"
    cache = tmp_path / "data"
    out.mkdir()
    cache.mkdir()
    for name, text in _BUILT.items():
        _write(out / name, text)
    (out / "subdir").mkdir()  # not a file: must not be listed
    _write(cache / "go-basic.obo", _OBO)
    return str(out), str(cache)


@pytest.fixture
def live(monkeypatch):
    """Fake the live version lookups; record every (url, dest, refresh) call."""
    calls = []
    with open(os.path.join(FIXTURES, "info_kegg.txt")) as f:
        info = f.read()

    def fake(url, dest, cache_dir, refresh=False):
        calls.append((url, dest, refresh))
        body = info if url.endswith("/info/kegg") else _REGISTRY
        path = os.path.join(cache_dir, dest)
        _write(path, body)
        return path

    monkeypatch.setattr(am, "cached_download", fake)
    return calls


def _by_file(manifest):
    return {f["file"]: f for f in manifest["files"]}


def test_lists_exactly_the_built_files_with_their_sha256(build):
    out, cache = build
    m = am.build_annotation_manifest(out, cache, gene_annotation_source="gaf")
    assert [f["file"] for f in m["files"]] == sorted(_BUILT)
    for f in m["files"]:
        assert f["sha256"] == sha256_file(os.path.join(out, f["file"]))


def test_writing_adds_only_the_manifest_and_never_lists_itself(build):
    out, cache = build
    before = {n: sha256_file(os.path.join(out, n)) for n in _BUILT}
    am.write_annotation_manifest(out, cache, gene_annotation_source="gaf")
    am.write_annotation_manifest(out, cache, gene_annotation_source="gaf")  # rerun
    assert set(os.listdir(out)) == set(_BUILT) | {"subdir", "manifest.json"}
    assert {n: sha256_file(os.path.join(out, n)) for n in _BUILT} == before
    with open(os.path.join(out, "manifest.json")) as f:
        m = json.load(f)
    assert "manifest.json" not in _by_file(m)


def test_go_ontology_and_gene_annotation_are_separate_records(build):
    out, cache = build
    m = am.build_annotation_manifest(out, cache, gene_annotation_source="gaf")
    assert m["go_ontology"] == {
        "url": am.GO_OBO_URL,
        "data_version": "releases/2026-01-01",
        "sha256": sha256_file(os.path.join(cache, "go-basic.obo")),
    }
    # GAF details stay in the companion file; nothing is copied or guessed.
    assert m["gene_annotation"] == {"source": "gaf", "details_file": "GO_provenance.json"}
    files = _by_file(m)
    assert files["GO2gene_BP.tab"]["sources"] == ["go_ontology", "gene_annotation"]
    assert files["GO_BP.gmt"]["sources"] == ["go_ontology", "gene_annotation"]
    assert files["GO2name_BP.tab"]["sources"] == ["go_ontology"]
    assert files["Annotation.tab"]["sources"] == ["gene_annotation"]
    assert files["KEGG_pathway2gene.tab"]["sources"] == ["kegg"]
    assert files["GO_provenance.json"]["sources"] == []


def test_kegg_records_the_per_database_dates_only(build, live):
    out, cache = build
    m = am.build_annotation_manifest(out, cache, query_versions=True)
    expected = parse_kegg_info_dates(os.path.join(FIXTURES, "info_kegg.txt"))
    assert m["kegg"] == {"kegg_db_dates": expected, "details_file": None}


def test_without_query_versions_live_versions_are_null_and_nothing_is_fetched(build, live):
    out, cache = build
    m = am.build_annotation_manifest(out, cache, gene_annotation_source="ensembl")
    assert m["kegg"]["kegg_db_dates"] is None
    assert m["gene_annotation"]["registry_display_name"] is None
    assert m["gene_annotation"]["registry_host"] == "https://www.ensembl.org"
    assert live == []
    # The exact, build-used GO ontology version is still recorded.
    assert m["go_ontology"]["data_version"] == "releases/2026-01-01"


def test_live_versions_are_refetched_into_manifest_only_cache_files(build, live):
    out, cache = build
    # A builder cache from an earlier run: must be neither read nor touched.
    stale = os.path.join(cache, "kegg_info_kegg.txt")
    _write(stale, "kegg\n\tpathway        1  2001/01/01\n")
    _write(os.path.join(cache, "manifest_biomart_registry.xml"), "<old other host/>")
    m = am.build_annotation_manifest(out, cache, gene_annotation_source="ensembl",
                                     ensembl_host="https://plants.ensembl.org/",
                                     ensembl_mart="plants_mart", query_versions=True)
    assert m["kegg"]["kegg_db_dates"]["pathway"] == "2026-07-03"
    assert m["gene_annotation"]["registry_display_name"] == "Ensembl Plants Genes 60"
    assert sorted(live) == sorted([
        ("https://rest.kegg.jp/info/kegg", "manifest_kegg_info.txt", True),
        ("https://plants.ensembl.org/biomart/martservice?type=registry",
         "manifest_biomart_registry.xml", True),
    ])
    with open(stale) as f:
        assert "2001/01/01" in f.read()


def test_kegg_details_file_points_at_kegg_provenance(build):
    out, cache = build
    _write(os.path.join(out, "KEGG_provenance.json"), "{}\n")
    m = am.build_annotation_manifest(out, cache)
    assert m["kegg"]["details_file"] == "KEGG_provenance.json"


@pytest.mark.parametrize("mart,expected", [
    ("ensembl", "Ensembl Genes 113"),
    ("plants_mart", "Ensembl Plants Genes 60"),
    ("plants", "Ensembl Plants Genes 60"),
    ("metazoa_mart", None),
])
def test_parse_biomart_registry(mart, expected):
    assert am.parse_biomart_registry(_REGISTRY, mart) == expected


def test_ensembl_records_host_and_registry_name(build, live):
    out, cache = build
    m = am.build_annotation_manifest(out, cache, gene_annotation_source="ensembl",
                                     query_versions=True)
    assert m["gene_annotation"]["registry_host"] == "https://www.ensembl.org"
    assert m["gene_annotation"]["registry_display_name"] == "Ensembl Genes 113"


def test_unavailable_version_queries_are_null_and_do_not_fail(build, monkeypatch):
    out, cache = build

    def boom(*args, **kwargs):
        raise RuntimeError("network down")

    monkeypatch.setattr(am, "cached_download", boom)
    path = am.write_annotation_manifest(out, cache, gene_annotation_source="ensembl",
                                        query_versions=True)
    with open(path) as f:
        m = json.load(f)
    assert m["kegg"]["kegg_db_dates"] is None
    assert m["gene_annotation"]["registry_display_name"] is None
    assert len(m["files"]) == len(_BUILT)


def test_unknown_versions_are_null(build):
    out, cache = build
    os.remove(os.path.join(cache, "go-basic.obo"))
    m = am.build_annotation_manifest(out, cache, gene_annotation_source="uniprot")
    assert m["go_ontology"] is None
    assert m["gene_annotation"] == {"source": "uniprot", "details_file": "GO_provenance.json"}

    _write(os.path.join(cache, "go-basic.obo"), "format-version: 1.2\n\n[Term]\n")
    assert am.build_annotation_manifest(out, cache)["go_ontology"]["data_version"] is None
    # No source given -> no gene_annotation record, rather than a guess.
    assert am.build_annotation_manifest(out, cache)["gene_annotation"] is None


def test_blocks_are_null_when_no_file_uses_them(tmp_path, build, live):
    _, cache = build
    out = tmp_path / "kegg_only"
    out.mkdir()
    _write(out / "KEGG_pathway2gene.tab", _BUILT["KEGG_pathway2gene.tab"])
    m = am.build_annotation_manifest(str(out), cache, gene_annotation_source="trinotate",
                                     query_versions=True)
    assert m["go_ontology"] is None and m["gene_annotation"] is None
    assert m["kegg"]["kegg_db_dates"]


def test_kegg_files_derive_from_ensembl_only_when_the_build_rekeyed_through_it(build):
    out, cache = build
    kegg = ["KEGG_pathway2gene.tab", "KEGG_pathway2name.tab"]

    rekeyed = _by_file(am.build_annotation_manifest(
        out, cache, gene_annotation_source="ensembl", kegg_ids_from_gene_annotation=True))
    for name in kegg:
        assert rekeyed[name]["sources"] == ["kegg", "gene_annotation"]
    # Non-KEGG files keep their own sources.
    assert rekeyed["GO2name_BP.tab"]["sources"] == ["go_ontology"]
    assert rekeyed["Annotation.tab"]["sources"] == ["gene_annotation"]

    # Same file names, other build paths (UniProt via KEGG /conv, no flag; or
    # no gene annotation at all): KEGG files derive from KEGG only.
    for kwargs in ({"gene_annotation_source": "uniprot"},
                   {"gene_annotation_source": None, "kegg_ids_from_gene_annotation": True}):
        files = _by_file(am.build_annotation_manifest(out, cache, **kwargs))
        for name in kegg:
            assert files[name]["sources"] == ["kegg"]


def test_workflows_record_inputs_as_submitted():
    import yaml
    for wf, expected in (("build-rna-annotation.yml",
                          ('"id_label=$ID_LABEL"', '"date=$IN_DATE"')),
                         ("build-model-annotation.yml", ('"date=$IN_DATE"',))):
        with open(os.path.join(ROOT, ".github", "workflows", wf)) as f:
            steps = yaml.safe_load(f)["jobs"]["build"]["steps"]
        run = next(s["run"] for s in steps if "run_annotation_manifest.py" in s.get("run", ""))
        for text in expected + ("--query-source-versions",):
            assert text in run, (wf, text)
        assert '"id_label=$label"' not in run


def test_inputs_and_release_are_recorded_verbatim(build):
    out, cache = build
    inputs = {"workflow": "build-rna-annotation", "project": "demo", "kegg_org": ""}
    m = am.build_annotation_manifest(out, cache, release_tag="rna_annot_demo_20260101",
                                     run_url="https://example.org/run/1", inputs=inputs)
    assert m["inputs"] == inputs
    assert m["release_tag"] == "rna_annot_demo_20260101"
    assert m["run_url"] == "https://example.org/run/1"
    assert set(m) == {"release_tag", "run_url", "builder_git_sha", "build_timestamp_utc",
                      "inputs", "go_ontology", "gene_annotation", "kegg", "files"}


def test_unknown_gene_annotation_source_is_rejected(build):
    out, cache = build
    with pytest.raises(ValueError):
        am.build_annotation_manifest(out, cache, gene_annotation_source="biomart")


def test_cli_writes_the_manifest(build):
    out, cache = build
    script = os.path.join(ROOT, "scripts", "run_annotation_manifest.py")
    subprocess.run([sys.executable, script, "--out", out, "--cache", cache,
                    "--release-tag", "t1", "--gene-annotation-source", "ensembl",
                    "--kegg-ids-from-gene-annotation",
                    "--input", "project=demo", "--input", "id_label="],
                   check=True, capture_output=True)
    with open(os.path.join(out, "manifest.json")) as f:
        m = json.load(f)
    assert m["inputs"] == {"project": "demo", "id_label": ""}
    assert m["gene_annotation"]["source"] == "ensembl"
    # No --query-source-versions: no live query, so no release name is claimed.
    assert m["gene_annotation"]["registry_display_name"] is None
    assert _by_file(m)["KEGG_pathway2gene.tab"]["sources"] == ["kegg", "gene_annotation"]

    bad = subprocess.run([sys.executable, script, "--out", out, "--cache", cache,
                          "--input", "no-equals-sign"], capture_output=True)
    assert bad.returncode != 0
