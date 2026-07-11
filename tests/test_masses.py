"""Tests for src/masses.py: formula validation + neutral monoisotopic mass.

The arithmetic is mass2chem's; these tests guard the validation layer that keeps
mass2chem from being fed formulas it silently mishandles.
"""

import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

pytest.importorskip("mass2chem", reason="mass2chem is a scoped optional dep")

from src.masses import is_computable_formula, neutral_mono_mass  # noqa: E402


# (formula, expected monoisotopic mass) -- independent reference values.
KNOWN = [
    ("C6H12O6", 180.063388),   # D-Glucose
    ("H2O", 18.010565),
    ("C3H4O3", 88.016044),     # pyruvate
    ("C10H16N5O13P3", 506.995747),  # ATP
    ("C21H27N7O14P2", 663.109123),  # NAD+
]


@pytest.mark.parametrize("formula,expected", KNOWN)
def test_known_masses_within_1mDa(formula, expected):
    m = neutral_mono_mass(formula)
    assert m is not None
    assert abs(m - expected) < 0.001  # within 1 mDa


@pytest.mark.parametrize("formula", [f for f, _ in KNOWN] + ["Fe4S4", "C16H30O4Zn"])
def test_valid_formulas_accepted(formula):
    assert is_computable_formula(formula) is True


@pytest.mark.parametrize("formula", [
    "(C6H10O5)n",          # polymer -> mass2chem silently returns monomer mass
    "C00031",              # a KEGG id, not a formula (leading-zero count)
    "C10H12N5O6PR2",       # generic R group
    "C6H10O5R",            # trailing R group
    "C2H12O-5",            # signed count
    "(C5H8NO4)n.C2H4O2",   # parenthesised polymer + hydrate
    "",
    None,
    "*",
])
def test_bad_formulas_rejected(formula):
    assert is_computable_formula(formula) is False
    assert neutral_mono_mass(formula) is None


def test_rounding_places():
    assert neutral_mono_mass("C6H12O6", decimals=3) == 180.063
