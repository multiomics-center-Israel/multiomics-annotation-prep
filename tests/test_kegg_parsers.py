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
    parse_org_pathways,
    parse_org_reaction_links,
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


def test_org_reaction_links():
    ids, gene2rxn = parse_org_reaction_links(
        os.path.join(FIX, "link_reaction_cre.txt"))
    assert "R00299" in ids
    assert ids.count("R00299") == 1  # de-duplicated
    assert "R00299" in gene2rxn["CHLRE_05g250000v5"]


def test_kegg_release_parsed():
    assert parse_kegg_release(os.path.join(FIX, "info_cre.txt")) == "KEGG Release 110.0"
