"""Tests for KEGG organism gene sets in src/prepare_kegg_org_sets.py.

Offline: the KEGG REST responses (link/pathway, list/pathway, info) are
pre-seeded in the cache dir under the names cached_download uses, so no
network is touched.
"""

import json
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.prepare_kegg_org_sets import (  # noqa: E402
    prepare_kegg_org_sets,
    remap_pathway_genes,
)

_LINKS = ("xyz:LOC_0001\tpath:xyz00010\n"
          "xyz:LOC_0002\tpath:xyz00010\n"
          "xyz:LOC_0002\tpath:xyz00020\n"
          "xyz:LOC_0003\tpath:xyz00020\n")
_NAMES = ("path:xyz00010\tGlycolysis - Test organism\n"
          "path:xyz00020\tTCA cycle - Test organism\n")
_INFO = "xyz              Test organism KEGG Genes Database\n" \
        "xyz              Release 999.0+/01-01, Jan 26\n"


def _seed(tmp_path):
    cache = tmp_path / "cache"
    cache.mkdir()
    (cache / "kegg_xyz_gene2path.txt").write_text(_LINKS)
    (cache / "kegg_xyz_pathway_names.txt").write_text(_NAMES)
    (cache / "kegg_xyz_info.txt").write_text(_INFO)
    return str(cache)


def _read_gmt(path):
    with open(path) as f:
        return {c[0]: (c[1], set(c[2:]))
                for c in (line.rstrip("\n").split("\t") for line in f)}


def test_remap_pathway_genes_fans_out_and_drops_unmapped():
    pairs = remap_pathway_genes(["p1", "p1", "p2"], ["A", "B", "C"],
                                {"A": {"A1", "A2"}, "C": {"C1"}})
    assert sorted(pairs) == [("p1", "A1"), ("p1", "A2"), ("p2", "C1")]


def test_remap_pathway_genes_errors_when_nothing_maps():
    with pytest.raises(RuntimeError, match="KEGG gene ids"):
        remap_pathway_genes(["p1"], ["A"], {"OTHER": {"X"}})


def test_kegg_org_sets_without_id_map_keeps_kegg_ids(tmp_path):
    cache = _seed(tmp_path)
    out = str(tmp_path / "out")
    prepare_kegg_org_sets("xyz", out, cache)
    gmt = _read_gmt(os.path.join(out, "KEGG_pathway.gmt"))
    assert gmt["xyz00010"] == ("Glycolysis - Test organism", {"LOC_0001", "LOC_0002"})
    assert gmt["xyz00020"][1] == {"LOC_0002", "LOC_0003"}
    with open(os.path.join(out, "KEGG_pathway2gene.tab")) as f:
        assert f.readline().rstrip("\n") == "v1\tindex"


def test_kegg_org_sets_with_id_map_and_provenance(tmp_path):
    cache = _seed(tmp_path)
    out = str(tmp_path / "out")
    id_map = tmp_path / "map.tsv"
    id_map.write_text("LOC_0001\tP1.1\nLOC_0002\tP2.1\nLOC_0002\tP2.2\n")
    prepare_kegg_org_sets("xyz", out, cache, id_map=str(id_map))

    gmt = _read_gmt(os.path.join(out, "KEGG_pathway.gmt"))
    assert gmt["xyz00010"][1] == {"P1.1", "P2.1", "P2.2"}
    assert gmt["xyz00020"][1] == {"P2.1", "P2.2"}  # LOC_0003 had no mapping
    with open(os.path.join(out, "KEGG_provenance.json")) as f:
        prov = json.load(f)
    assert prov["kegg_org"] == "xyz"
    assert "Release 999.0+/01-01, Jan 26" in " ".join(prov["kegg_release"])
    assert prov["id_map"]["source"] == "file"
    assert len(prov["links"]["sha256"]) == 64


@pytest.mark.parametrize("bad", ["GLA", "g1a", "giardia", "../x"])
def test_kegg_org_sets_rejects_bad_org_code(tmp_path, bad):
    with pytest.raises(ValueError):
        prepare_kegg_org_sets(bad, str(tmp_path), str(tmp_path))
