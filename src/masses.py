"""Neutral monoisotopic mass for the mummichog metabolic-model module.

The mummichog model needs each compound's *neutral* monoisotopic mass (mummichog
adds adducts itself). We do NOT reimplement mass calculation: the arithmetic is
delegated to :func:`mass2chem.formula.calculate_formula_mass`.

This module only decides whether a KEGG ``FORMULA`` string is a *concrete neutral
formula* that mass2chem can be trusted with. mass2chem's parser
(``([A-Z][a-z]*)(\\d*)``) silently mis-handles a few shapes that appear in KEGG:

* polymers / groups in parentheses, e.g. ``(C6H10O5)n`` -> monomer mass, no error
* repeated element symbols -> later occurrence overwrites the earlier count
* signed counts, e.g. ``C2H12O-5`` -> wrong dict, no error
* generic ``R`` / ``X`` groups -> ``KeyError`` on an unknown element

so we reject those up front and only pass through clean formulas.

These are optional dependencies scoped to this module; importing ``src.masses``
does not require them until a mass is actually computed, keeping the gene-set
modules dependency-free.
"""

import re

# Element symbol + optional count, matching mass2chem's own tokenizer so that a
# formula we accept is tokenized identically by mass2chem.
_TOKEN = re.compile(r"([A-Z][a-z]*)(\d*)")
_ALNUM = re.compile(r"[A-Za-z0-9]+")


def _atom_mass_dict():
    """Return mass2chem's element mass table (lazy optional import)."""
    try:
        from mass2chem.formula import atom_mass_dict
    except ImportError as exc:  # pragma: no cover - exercised via install hint
        raise ImportError(
            "mass2chem is required for the mummichog model module. Install the "
            "scoped optional deps:  pip install -r requirements-mummichog.txt"
        ) from exc
    return atom_mass_dict


def is_computable_formula(formula):
    """True iff *formula* is a concrete neutral formula mass2chem handles correctly.

    Rejects empty/None, anything with characters other than ``[A-Za-z0-9]``
    (parentheses, dots, +/- charges, ``*`` wildcards), leading-zero counts
    (which is what a KEGG id such as ``C00031`` looks like), repeated element
    symbols, and unknown elements (``R``/``X`` generic groups and glycan
    placeholders).
    """
    if not formula:
        return False
    f = formula.strip()
    if not _ALNUM.fullmatch(f):
        return False
    tokens = _TOKEN.findall(f)
    # Reconstruction must be lossless (no leftover characters the tokenizer drops).
    if "".join(el + n for el, n in tokens) != f:
        return False
    atom_mass = _atom_mass_dict()
    seen = set()
    for el, n in tokens:
        if el not in atom_mass:
            return False
        if el in seen:
            # duplicate symbol -> mass2chem overwrites the count; not trustworthy
            return False
        seen.add(el)
        if n and n.lstrip("0") != n:
            # leading-zero count, e.g. the "00031" in a KEGG id "C00031"
            return False
    return True


def neutral_mono_mass(formula, decimals=6):
    """Monoisotopic mass of a neutral *formula*, or ``None`` if not computable.

    The actual summation is done by mass2chem; this function is a thin, validated
    wrapper so callers never feed mass2chem a formula it would mishandle.
    """
    if not is_computable_formula(formula):
        return None
    try:
        from mass2chem.formula import calculate_formula_mass
    except ImportError as exc:  # pragma: no cover - exercised via install hint
        raise ImportError(
            "mass2chem is required for the mummichog model module. Install the "
            "scoped optional deps:  pip install -r requirements-mummichog.txt"
        ) from exc
    try:
        return round(calculate_formula_mass(formula.strip()), decimals)
    except Exception:  # noqa: BLE001 - any parse/lookup failure -> not computable
        return None
