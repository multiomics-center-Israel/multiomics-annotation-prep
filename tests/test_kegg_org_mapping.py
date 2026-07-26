"""Tests for the model-organism KEGG mapping in src/download_kegg_org.py.

These run offline: KEGG downloads are stubbed by writing the fixture-shaped
files into a tmp cache and monkeypatching ``cached_download`` to return them.
No network, no heavy deps.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src import download_kegg_org as dko  # noqa: E402
from src.download_kegg_org import (  # noqa: E402
    kegg_to_ensembl_map,
    prepare_kegg_by_org,
)


def _write(path, text):
    with open(path, "w") as f:
        f.write(text)


def _stub_downloads(monkeypatch, files):
    """Route cached_download(url, name, ...) to pre-written files by name."""
    def fake(url, name, cache_dir, refresh=False, **kw):
        return files[name]
    monkeypatch.setattr(dko, "cached_download", fake)


def test_kegg_to_ensembl_map_bridges_via_ncbi(monkeypatch, tmp_path):
    # KEGG gene -> NCBI gene id (from conv), NCBI -> Ensembl (from BioMart).
    monkeypatch.setattr(dko, "kegg_conv_map", lambda *a, **k: {
        "mmu:11298": "11298",   # has an Ensembl xref
        "mmu:11302": "11302",   # has an Ensembl xref
        "mmu:99999": "99999",   # NCBI id with no Ensembl xref -> dropped
    })
    ncbi_to_ensembl = {
        "11298": "ENSMUSG00000000001",
        "11302": "ENSMUSG00000000002",
    }
    m = kegg_to_ensembl_map("mmu", ncbi_to_ensembl, str(tmp_path))
    assert m == {
        "mmu:11298": "ENSMUSG00000000001",
        "mmu:11302": "ENSMUSG00000000002",
    }
    assert "mmu:99999" not in m  # unmapped NCBI id is not carried over


def test_prepare_kegg_by_org_ext_map_and_direct_fallback(monkeypatch, tmp_path):
    # Two genes: one resolved by the explicit bridge map, one only by the
    # direct KEGG-id == Ensembl-id fallback (the Arabidopsis-style case), and
    # one that neither resolves (dropped).
    gene2path = tmp_path / "g2p.txt"
    pwnames = tmp_path / "pwn.txt"
    _write(gene2path, "\t".join(["ath:AT1G01010", "path:ath00010"]) + "\n"
           + "\t".join(["ath:AT1G01020", "path:ath00010"]) + "\n"
           + "\t".join(["ath:AT9G99999", "path:ath00020"]) + "\n")
    _write(pwnames,
           "path:ath00010\tGlycolysis - Arabidopsis thaliana\n"
           "path:ath00020\tCitrate cycle - Arabidopsis thaliana\n")
    _stub_downloads(monkeypatch, {
        "kegg_ath_gene2path.txt": str(gene2path),
        "kegg_ath_pathway_names.txt": str(pwnames),
    })

    out = tmp_path / "out"
    out.mkdir()
    ext_id_map = {"ath:AT1G01010": "AT1G01010"}       # via ncbi bridge
    ext_universe = {"AT1G01010", "AT1G01020"}          # direct fallback set

    prepare_kegg_by_org("ath", str(out), str(tmp_path), id_source="ensembl",
                        ext_id_map=ext_id_map, ext_id_universe=ext_universe)

    with open(out / "KEGG_pathway2gene.tab") as f:
        body = f.read()
    assert "AT1G01010" in body            # from ext_id_map
    assert "AT1G01020" in body            # from direct fallback
    assert "AT9G99999" not in body        # neither -> dropped
    assert "ath:" not in body             # prefixes stripped everywhere

    with open(out / "KEGG_pathway2name.tab") as f:
        names = f.read()
    assert "ath00010" in names
    assert "Glycolysis" in names


def test_prepare_kegg_by_org_rejects_ensembl_without_bridge(monkeypatch, tmp_path):
    # Guard: id_source='ensembl' with no ext map must fail fast, not hit the
    # invalid conv/ensembl/<org> endpoint.
    gene2path = tmp_path / "g2p.txt"
    pwnames = tmp_path / "pwn.txt"
    _write(gene2path, "ath:AT1G01010\tpath:ath00010\n")
    _write(pwnames, "path:ath00010\tGlycolysis\n")
    _stub_downloads(monkeypatch, {
        "kegg_ath_gene2path.txt": str(gene2path),
        "kegg_ath_pathway_names.txt": str(pwnames),
    })
    out = tmp_path / "out"
    out.mkdir()
    try:
        prepare_kegg_by_org("ath", str(out), str(tmp_path), id_source="ensembl")
    except ValueError as exc:
        assert "ensembl" in str(exc).lower()
    else:
        raise AssertionError("expected ValueError for id_source='ensembl' "
                             "without a bridge map")


def test_prepare_kegg_by_org_conv_uniprot(monkeypatch, tmp_path):
    # id_source='uniprot' uses KEGG /conv/uniprot directly (a valid conv DB).
    gene2path = tmp_path / "g2p.txt"
    pwnames = tmp_path / "pwn.txt"
    conv = tmp_path / "conv.txt"
    _write(gene2path, "hsa:5091\tpath:hsa00010\n" "hsa:5162\tpath:hsa00020\n")
    _write(pwnames, "path:hsa00010\tGlycolysis\npath:hsa00020\tTCA\n")
    _write(conv, "hsa:5091\tup:P35557\nhsa:5162\tup:P11177\n")
    _stub_downloads(monkeypatch, {
        "kegg_hsa_gene2path.txt": str(gene2path),
        "kegg_hsa_pathway_names.txt": str(pwnames),
        "kegg_hsa_to_uniprot.txt": str(conv),
    })
    out = tmp_path / "out"
    out.mkdir()
    prepare_kegg_by_org("hsa", str(out), str(tmp_path), id_source="uniprot")
    with open(out / "KEGG_pathway2gene.tab") as f:
        body = f.read()
    assert "P35557" in body and "P11177" in body
    assert "hsa:" not in body and "up:" not in body
