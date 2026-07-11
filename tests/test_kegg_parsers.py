"""Tests for the KEGG compound/reaction/pathway parsers in src/download_kegg.py.

These run against small saved fixtures (no network) and do not need the optional
mass2chem / metDataModel deps.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.download_kegg import (  # noqa: E402
    iter_kegg_records,
    parse_compound_record,
    parse_kegg_release,
    parse_ko_reaction_links,
    parse_org_ko_links,
    parse_org_pathways,
    parse_reaction_record,
)

FIX = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures")


def _records(name):
    with open(os.path.join(FIX, name)) as f:
        return list(iter_kegg_records(f.read()))


def test_compound_records_parsed():
    cpds = {c["id"]: c for c in (parse_compound_record(r)
                                 for r in _records("kegg_compounds.txt"))}
    assert cpds["C00031"]["name"] == "D-Glucose"      # first name, ';' stripped
    assert cpds["C00031"]["formula"] == "C6H12O6"
    assert cpds["C00031"]["exact_mass"] == 180.0634
    # polymer compound: formula present but no exact mass
    assert cpds["C00369"]["formula"] == "(C6H10O5)n"
    assert cpds["C00369"]["exact_mass"] is None


def test_reaction_equation_splits_into_substrates_and_products():
    rxns = {r["id"]: r for r in (parse_reaction_record(r)
                                 for r in _records("kegg_reactions.txt"))}
    r = rxns["R00299"]
    assert r["reactants"] == ["C00002", "C00031"]
    assert r["products"] == ["C00008", "C00668"]
    assert r["enzymes"] == ["2.7.1.1", "2.7.1.2", "2.7.1.63"]
    assert "00010" in r["ref_pathways"]


def test_reaction_ec_and_pathway_extraction():
    rxns = {r["id"]: r for r in (parse_reaction_record(r)
                                 for r in _records("kegg_reactions.txt"))}
    assert rxns["R00771"]["enzymes"] == ["5.3.1.9"]
    assert rxns["R00771"]["reactants"] == ["C00668"]
    assert rxns["R00771"]["products"] == ["C00085"]


def test_org_pathways_strip_org_suffix():
    pw = parse_org_pathways(os.path.join(FIX, "cre_pathways.txt"))
    assert pw["cre00010"] == "Glycolysis / Gluconeogenesis"
    assert "cre00020" in pw


def test_org_ko_links():
    kos, gene2ko = parse_org_ko_links(os.path.join(FIX, "link_ko_cre.txt"))
    assert "K00844" in kos
    assert kos.count("K00844") == 1  # de-duplicated across genes
    assert gene2ko["CHLRE_01g000050v5"] == ["K00844"]


def test_ko_to_reaction_step():
    # Guards the KO->reaction resolution (organism genes are NOT linked directly
    # to reactions; they route through KO).
    ko2rxn = parse_ko_reaction_links(os.path.join(FIX, "link_reaction_ko.txt"))
    assert ko2rxn["K00844"] == ["R00299"]
    assert ko2rxn["K01810"] == ["R00771", "R01786"]
    # intersect a KO set with the KEGG-wide map -> only in-set KOs contribute
    org_kos = ["K00844", "K01810", "K00873", "K00688", "K99999"]
    reactions = sorted({r for ko in org_kos for r in ko2rxn.get(ko, [])})
    assert reactions == ["R00200", "R00299", "R00771", "R01786", "R02110"]
    assert "R99999" not in reactions  # its KO (K00001) is not in the organism
    assert ko2rxn.get("K99999") is None  # coverage gap


def test_kegg_release_parsed():
    # canonical release comes from info/kegg
    assert parse_kegg_release(os.path.join(FIX, "info_kegg.txt")) == "KEGG Release 116.0"
