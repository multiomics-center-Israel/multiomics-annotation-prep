"""Tests for the GAF GO path in src/prepare_gaf.py.

Offline: a tiny synthetic GAF, id maps and go-basic.obo are written to tmp_path;
the OBO is pre-seeded in the cache dir so cached_download never hits the network.
"""

import gzip
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.prepare_gaf import (  # noqa: E402
    parse_gaf,
    parse_id_map,
    prepare_gaf,
    remap_genes,
)

_OBO = """format-version: 1.2

[Term]
id: GO:0000001
name: root process
namespace: biological_process

[Term]
id: GO:0000002
name: leaf process
namespace: biological_process
is_a: GO:0000001 ! root process

[Term]
id: GO:0000010
name: some function
namespace: molecular_function
"""

# GENE_B's only row is NOT-qualified, so GENE_B must not appear at all.
_GAF = "\n".join([
    "!gaf-version: 2.2",
    "\t".join(["TestDB", "GENE_A", "symA", "involved_in", "GO:0000002", "PMID:1",
               "IEA", "", "P", "protein A", "", "protein", "taxon:1", "20240101",
               "TestDB"]),
    "\t".join(["TestDB", "GENE_A", "symA", "enables", "GO:0000010", "PMID:1",
               "IEA", "", "F", "protein A", "", "protein", "taxon:1", "20240101",
               "TestDB"]),
    "\t".join(["TestDB", "GENE_B", "symB", "NOT|involved_in", "GO:0000002", "PMID:1",
               "IEA", "", "P", "protein B", "", "protein", "taxon:1", "20240101",
               "TestDB"]),
    "\t".join(["TestDB", "GENE_C", "", "involved_in", "GO:0000001", "PMID:1",
               "IEA", "", "P", "protein C", "", "protein", "taxon:1", "20240101",
               "TestDB"]),
]) + "\n"

_FEATURE_TABLE = "\n".join([
    "# feature\tclass\tassembly\tassembly_unit\tseq_type\tchromosome\tgenomic_accession"
    "\tstart\tend\tstrand\tproduct_accession\tnon-redundant_refseq\trelated_accession"
    "\tname\tsymbol\tGeneID\tlocus_tag\tfeature_interval_length\tproduct_length"
    "\tattributes",
    "gene\tprotein_coding\tASM1\tPrimary\tchromosome\t1\tCHR1\t1\t90\t+\t\t\t\t\t\t1"
    "\tGENE_A\t90\t\t",
    "CDS\twith_protein\tASM1\tPrimary\tchromosome\t1\tCHR1\t1\t90\t+\tPROT1.1\t\t\tp A"
    "\t\t1\tGENE_A\t90\t29\t",
    "CDS\twith_protein\tASM1\tPrimary\tchromosome\t1\tCHR1\t1\t90\t+\tPROT1.2\t\t\tp A"
    "\t\t1\tGENE_A\t90\t29\t",
]) + "\n"


def _write(path, text, gz=False):
    opener = gzip.open if gz else open
    with opener(path, "wt") as f:
        f.write(text)
    return str(path)


def _cache(tmp_path):
    cache = tmp_path / "cache"
    cache.mkdir()
    _write(cache / "go-basic.obo", _OBO)
    return str(cache)


def _read_gmt(path):
    out = {}
    with open(path) as f:
        for line in f:
            cols = line.rstrip("\n").split("\t")
            out[cols[0]] = (cols[1], set(cols[2:]))
    return out


@pytest.mark.parametrize("gz", [False, True])
def test_parse_gaf_drops_not_and_reads_info(tmp_path, gz):
    path = _write(tmp_path / "x.gaf", _GAF, gz=gz)
    gene2go, info = parse_gaf(path)
    assert gene2go == {"GENE_A": {"GO:0000002", "GO:0000010"},
                       "GENE_C": {"GO:0000001"}}
    assert info["GENE_A"] == ("symA", "protein A")


def test_parse_id_map_feature_table_uses_cds_rows(tmp_path):
    path = _write(tmp_path / "ft.txt.gz", _FEATURE_TABLE, gz=True)
    assert parse_id_map(path) == {"GENE_A": {"PROT1.1", "PROT1.2"}}


def test_parse_id_map_two_column(tmp_path):
    path = _write(tmp_path / "map.tsv", "# from\tto\nGENE_A\tP1\nGENE_C\tP3\n\n")
    assert parse_id_map(path) == {"GENE_A": {"P1"}, "GENE_C": {"P3"}}


def test_remap_genes_fans_out_and_drops_unmapped():
    gene2go = {"GENE_A": {"GO:1"}, "GENE_C": {"GO:2"}}
    info = {"GENE_A": ("a", "A"), "GENE_C": ("c", "C")}
    new_go, new_info = remap_genes(gene2go, info, {"GENE_A": {"P1", "P2"}})
    assert new_go == {"P1": {"GO:1"}, "P2": {"GO:1"}}
    assert new_info["P2"] == ("a", "A")


def test_remap_genes_errors_when_nothing_maps():
    with pytest.raises(RuntimeError):
        remap_genes({"GENE_A": {"GO:1"}}, {}, {"OTHER": {"P1"}})


def test_prepare_gaf_writes_expanded_gmt(tmp_path):
    cache = _cache(tmp_path)
    out = str(tmp_path / "out")
    gaf = _write(tmp_path / "x.gaf", _GAF)
    prepare_gaf(gaf, out, cache)

    bp = _read_gmt(os.path.join(out, "GO_BP.gmt"))
    # GENE_A's leaf annotation propagates to the root term
    assert bp["GO:0000002"] == ("leaf process", {"GENE_A"})
    assert bp["GO:0000001"] == ("root process", {"GENE_A", "GENE_C"})
    mf = _read_gmt(os.path.join(out, "GO_MF.gmt"))
    assert mf == {"GO:0000010": ("some function", {"GENE_A"})}
    for name in ("GO2gene_BP.tab", "GO2name_BP.tab", "Annotation.tab"):
        assert os.path.exists(os.path.join(out, name))
    assert not os.path.exists(os.path.join(out, "GO_CC.gmt"))


def test_prepare_gaf_with_id_map(tmp_path):
    cache = _cache(tmp_path)
    out = str(tmp_path / "out")
    gaf = _write(tmp_path / "x.gaf", _GAF)
    ft = _write(tmp_path / "ft.txt.gz", _FEATURE_TABLE, gz=True)
    prepare_gaf(gaf, out, cache, id_map=ft)

    bp = _read_gmt(os.path.join(out, "GO_BP.gmt"))
    assert bp["GO:0000001"][1] == {"PROT1.1", "PROT1.2"}
    with open(os.path.join(out, "Annotation.tab")) as f:
        rows = [line.rstrip("\n").split("\t") for line in f]
    assert rows[0] == ["Gene", "gene_name", "protein_name"]
    assert [r[0] for r in rows[1:]] == ["PROT1.1", "PROT1.2"]


def test_prepare_gaf_rejects_a_non_gaf_file(tmp_path):
    cache = _cache(tmp_path)
    table = _write(tmp_path / "trinotate_go.txt", "GENE_A\tGO:0000002\n")
    with pytest.raises(RuntimeError, match="GAF"):
        prepare_gaf(table, str(tmp_path / "out"), cache)
