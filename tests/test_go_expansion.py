"""Tests for GO OBO parsing + hierarchy expansion in src/download_go.py.

Offline: parses a tiny synthetic go-basic.obo written to tmp_path. Guards that
the ancestor closure follows BOTH is_a and part_of (matching GO.db ANCESTOR /
clusterProfiler::buildGOmap), while excluding regulates and has_part.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.download_go import (  # noqa: E402
    build_ancestor_cache,
    get_ancestors,
    parse_obo,
)

# root -is_a- mid -is_a- leaf  (BP chain)
# whole <-part_of- part        (CC part_of edge)
# regulator -regulates-> root  (must NOT become an ancestor)
# haspart_term -has_part-> part (inverse edge; must NOT become an ancestor)
_OBO = """format-version: 1.2

[Term]
id: GO:0000001
name: root process
namespace: biological_process

[Term]
id: GO:0000002
name: mid process
namespace: biological_process
is_a: GO:0000001 ! root process

[Term]
id: GO:0000003
name: leaf process
namespace: biological_process
is_a: GO:0000002 ! mid process

[Term]
id: GO:0000010
name: whole structure
namespace: cellular_component

[Term]
id: GO:0000011
name: part structure
namespace: cellular_component
is_a: GO:0000012 ! generic part
relationship: part_of GO:0000010 ! whole structure

[Term]
id: GO:0000012
name: generic part
namespace: cellular_component

[Term]
id: GO:0000020
name: regulator process
namespace: biological_process
relationship: regulates GO:0000001 ! root process

[Term]
id: GO:0000030
name: assembly holder
namespace: cellular_component
relationship: has_part GO:0000011 ! part structure

[Term]
id: GO:0000099
name: obsolete thing
namespace: biological_process
is_a: GO:0000001 ! root process
is_obsolete: true
"""


def _parse(tmp_path):
    obo = tmp_path / "go-basic.obo"
    obo.write_text(_OBO)
    return parse_obo(str(obo))


def test_terms_and_namespaces(tmp_path):
    terms, _ = _parse(tmp_path)
    assert terms["GO:0000001"] == {"name": "root process", "namespace": "BP"}
    assert terms["GO:0000010"]["namespace"] == "CC"
    # obsolete terms are dropped
    assert "GO:0000099" not in terms


def test_is_a_chain_expands_transitively(tmp_path):
    _, parents = _parse(tmp_path)
    anc = get_ancestors("GO:0000003", parents)
    assert anc == {"GO:0000002", "GO:0000001"}


def test_part_of_is_followed(tmp_path):
    _, parents = _parse(tmp_path)
    # GO:0000011 has is_a GO:0000012 AND part_of GO:0000010 -> both are ancestors
    anc = get_ancestors("GO:0000011", parents)
    assert anc == {"GO:0000012", "GO:0000010"}


def test_regulates_is_not_an_ancestor(tmp_path):
    _, parents = _parse(tmp_path)
    # regulator -regulates-> root must not propagate
    assert get_ancestors("GO:0000020", parents) == set()
    assert "GO:0000020" not in parents


def test_has_part_is_not_an_ancestor(tmp_path):
    _, parents = _parse(tmp_path)
    # has_part is the inverse of part_of and points to a child, not a parent
    assert get_ancestors("GO:0000030", parents) == set()
    assert "GO:0000030" not in parents


def test_ancestor_cache_matches_direct(tmp_path):
    _, parents = _parse(tmp_path)
    cache = build_ancestor_cache(parents)
    assert get_ancestors("GO:0000011", parents, cache) == {"GO:0000012", "GO:0000010"}
    assert get_ancestors("GO:0000003", parents, cache) == {"GO:0000002", "GO:0000001"}
