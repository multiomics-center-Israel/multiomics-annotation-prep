"""Tests for the KEGG reaction-pair reference builder.

Parsing/pair-generation tests run against small saved fixtures (no network).
The orchestration tests monkeypatch the KEGG fetch layer so the fatal
fetch_failed path and the atomic-finalize happy path are exercised offline.
"""

import json
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import src.prepare_kegg_reaction_pairs as mod  # noqa: E402
from src.download_kegg import iter_kegg_records, parse_reaction_record  # noqa: E402
from src.prepare_kegg_reaction_pairs import build_reaction_pairs  # noqa: E402
from src.utils import sha256_file  # noqa: E402

FIX = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures")


def _reactions(name):
    with open(os.path.join(FIX, name)) as f:
        return [parse_reaction_record(r) for r in iter_kegg_records(f.read())]


# --- parsing: equation + arrow retained (corr. #1) -------------------------

def test_parse_reaction_record_retains_equation_and_arrow():
    rxns = {r["id"]: r for r in _reactions("kegg_reactions.txt")}
    r = rxns["R00299"]
    assert r["equation"] == "C00002 + C00031 <=> C00008 + C00668"
    assert r["equation_arrow"] == "<=>"
    # existing keys still present (backward compatible)
    assert r["reactants"] == ["C00002", "C00031"]
    assert r["products"] == ["C00008", "C00668"]


# --- cross-side Cartesian product with explicit expected rows --------------

def test_build_reaction_pairs_cross_side_expected_rows():
    rows, _ = build_reaction_pairs(_reactions("kegg_reactions.txt"))
    got = sorted((s, p) for rid, s, p, _e, _a in rows if rid == "R00299")
    assert got == [("C00002", "C00008"), ("C00002", "C00668"),
                   ("C00031", "C00008"), ("C00031", "C00668")]
    # single-pair reaction
    got_771 = [(s, p) for rid, s, p, _e, _a in rows if rid == "R00771"]
    assert got_771 == [("C00668", "C00085")]


def test_co_substrate_and_co_product_pairs_excluded():
    rows, _ = build_reaction_pairs(_reactions("kegg_reactions.txt"))
    pairs = {(s, p) for rid, s, p, _e, _a in rows if rid == "R00299"}
    # co-substrates (both left) never connected, in either direction
    assert ("C00002", "C00031") not in pairs
    assert ("C00031", "C00002") not in pairs
    # co-products (both right) never connected either
    assert ("C00008", "C00668") not in pairs
    assert ("C00668", "C00008") not in pairs


def test_duplicate_pair_from_two_reactions_preserved():
    rows, _ = build_reaction_pairs(_reactions("reaction_pairs_cases.txt"))
    dup = sorted(rid for rid, s, p, _e, _a in rows
                 if (s, p) == ("C10001", "C10002"))
    assert dup == ["R99005", "R99006"]  # kept distinct, not deduped to one


# --- precise exclusion reason codes (corr. #8) -----------------------------

def test_exclusion_reason_codes():
    _, excluded = build_reaction_pairs(_reactions("reaction_pairs_cases.txt"))
    reasons = {rid: reason for rid, reason, _ctx in excluded}
    assert reasons["R99001"] == "only_self_pairs"
    assert reasons["R99002"] == "missing_equation"
    assert reasons["R99003"] == "empty_left_compound_side"
    assert reasons["R99004"] == "unsupported_equation_arrow"
    # the two duplicate-pair reactions are NOT excluded
    assert "R99005" not in reasons and "R99006" not in reasons


# --- deterministic ordering (corr. #4) -------------------------------------

def test_deterministic_ordering():
    rxns = _reactions("kegg_reactions.txt")
    assert build_reaction_pairs(rxns)[0] == build_reaction_pairs(rxns)[0]


# --- orchestration: fetch_failed is fatal, writes no asset (corr. #3) ------

_REC = {"id": "R00001", "equation": "C00001 <=> C00002", "equation_arrow": "<=>",
        "reactants": ["C00001"], "products": ["C00002"]}


def _patch_fetch(monkeypatch, parsed, missing):
    monkeypatch.setattr(mod, "download_reaction_list", lambda c, refresh=False: "x")
    monkeypatch.setattr(mod, "parse_reaction_list",
                        lambda p: sorted(set(parsed) | set(missing)))
    monkeypatch.setattr(mod, "fetch_and_reconcile_reactions",
                        lambda ids, cache, **k: (dict(parsed), list(missing)))
    monkeypatch.setattr(mod, "_kegg_provenance", lambda cache, refresh: {
        "source_snapshot": "KEGG snapshot 2026-07-07",
        "reaction_database_date": "2026-07-07",
        "retrieved_at": "2026-07-13T00:00:00Z", "kegg_info_raw": None})


def test_fetch_failed_is_fatal_and_writes_no_asset(tmp_path, monkeypatch):
    _patch_fetch(monkeypatch, {"R00001": _REC}, ["R00002"])
    out = tmp_path / "out"
    out.mkdir()
    with pytest.raises(RuntimeError, match="fetch_failed"):
        mod.prepare_kegg_reaction_pairs(str(out), str(tmp_path / "cache"))
    names = [p.name for p in out.iterdir()]
    assert not any(n.endswith(".tsv.gz") for n in names)
    assert not any(n.endswith(".manifest.json") for n in names)
    assert any(n.endswith(".excluded.tsv") for n in names)  # failure report only


def test_allow_incomplete_produces_no_asset(tmp_path, monkeypatch):
    _patch_fetch(monkeypatch, {"R00001": _REC}, ["R00002"])
    out = tmp_path / "out"
    out.mkdir()
    data, manifest, excl = mod.prepare_kegg_reaction_pairs(
        str(out), str(tmp_path / "cache"), allow_incomplete=True)
    assert data is None and manifest is None
    assert excl.endswith(".excluded.tsv") and os.path.exists(excl)


# --- orchestration: complete build finalizes an asset atomically -----------

def test_prepare_writes_checksummed_asset_when_complete(tmp_path, monkeypatch):
    recs = {r["id"]: r for r in _reactions("kegg_reactions.txt")}
    _patch_fetch(monkeypatch, recs, [])
    out = tmp_path / "out"
    out.mkdir()
    data, manifest, excl = mod.prepare_kegg_reaction_pairs(
        str(out), str(tmp_path / "cache"))

    assert os.path.exists(data) and data.endswith(".tsv.gz")
    m = json.load(open(manifest))
    assert m["schema_version"] == 1
    assert m["pair_definition_method"] == "equation_side_cartesian_product"
    assert m["counts"]["n_rows"] > 0
    assert m["sha256"] == sha256_file(data)   # manifest checksum matches file
    # no staging dir left behind
    assert not any(p.name.startswith(".staging") for p in out.iterdir())
