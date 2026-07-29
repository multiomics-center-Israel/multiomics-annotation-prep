"""Tests for the direct BioMart REST layer in src/prepare_ensembl.py.

Offline: cached_download is stubbed to return canned martservice responses, so
no network and no pybiomart/pandas are needed.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src import prepare_ensembl as pe  # noqa: E402


def _stub(monkeypatch, content):
    def _cd(url, dest, cache_dir, refresh=False):
        path = os.path.join(cache_dir, dest)
        with open(path, "w", encoding="utf-8") as f:
            f.write(content)
        return path
    monkeypatch.setattr(pe, "cached_download", _cd)


def test_virtual_schema_mapping():
    assert pe._virtual_schema_for("ensembl") == "default"
    assert pe._virtual_schema_for(None) == "default"
    assert pe._virtual_schema_for("default") == "default"
    assert pe._virtual_schema_for("plants") == "plants_mart"
    assert pe._virtual_schema_for("plants_mart") == "plants_mart"
    assert pe._virtual_schema_for("protists_mart") == "protists_mart"


def test_biomart_query_parses_tsv_and_pads(monkeypatch, tmp_path):
    _stub(monkeypatch,
          "Gene\tName\tBiotype\tDesc\n"
          "G1\tnameA\tprotein_coding\tdesc one\n"
          "G2\t\tncRNA\t\n")            # trailing empty fields
    rows = pe._biomart_query("https://x", "ds", ["a", "b", "c", "d"],
                             "default", str(tmp_path), "t.tsv")
    assert rows[0] == ["Gene", "Name", "Biotype", "Desc"]
    assert rows[1] == ["G1", "nameA", "protein_coding", "desc one"]
    assert rows[2] == ["G2", "", "ncRNA", ""]   # padded to 4 columns


def test_biomart_query_rejects_html_and_clears_cache(monkeypatch, tmp_path):
    _stub(monkeypatch, "<!DOCTYPE html><html><body>maintenance</body></html>")
    try:
        pe._biomart_query("https://x", "ds", ["a", "b"], "default",
                          str(tmp_path), "bad.tsv")
    except RuntimeError as exc:
        assert "non-TSV" in str(exc)
        assert not (tmp_path / "bad.tsv").exists()  # poisoned cache removed
    else:
        raise AssertionError("expected RuntimeError for an HTML response")


def test_biomart_query_rejects_query_error(monkeypatch, tmp_path):
    _stub(monkeypatch, "Query ERROR: Attribute foo not found\n")
    try:
        pe._biomart_query("https://x", "ds", ["ensembl_gene_id", "foo"],
                          "default", str(tmp_path), "err.tsv")
    except RuntimeError:
        pass
    else:
        raise AssertionError("expected RuntimeError for a Query ERROR response")


def test_fetch_ncbi_xref_falls_through_attrs(monkeypatch, tmp_path):
    # First attribute name errors; second returns a mapping.
    calls = {"n": 0}

    def _cd(url, dest, cache_dir, refresh=False):
        path = os.path.join(cache_dir, dest)
        calls["n"] += 1
        if "entrezgene_id" in dest:
            body = "Query ERROR: Attribute entrezgene_id not found\n"
        else:  # entrezgene
            body = "ensembl_gene_id\tentrezgene\nENSG1\t111\nENSG2\t222\n"
        with open(path, "w", encoding="utf-8") as f:
            f.write(body)
        return path
    monkeypatch.setattr(pe, "cached_download", _cd)

    m = pe._fetch_ncbi_xref("https://x", "hsapiens_gene_ensembl", "default",
                            str(tmp_path))
    assert m == {"111": "ENSG1", "222": "ENSG2"}
