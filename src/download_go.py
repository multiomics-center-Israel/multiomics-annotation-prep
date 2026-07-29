"""GO term names, namespaces, and OBO parsing."""

from collections import defaultdict

from .utils import cached_download, log_msg

NS_MAP = {
    "biological_process": "BP",
    "molecular_function": "MF",
    "cellular_component": "CC",
}

GO_OBO_URL = "http://purl.obolibrary.org/obo/go/go-basic.obo"


def go_term_table(cache_dir, refresh=False):
    """Return dict: go_id -> {name, namespace} and the ancestor graph."""
    obo_path = cached_download(GO_OBO_URL, "go-basic.obo", cache_dir, refresh)
    terms, parents = parse_obo(obo_path)
    return terms, parents


def parse_obo(path):
    """Parse go-basic.obo -> (terms dict, parents dict).

    terms: go_id -> {name, namespace}
    parents: go_id -> set of parent go_ids

    Parents include both ``is_a`` and ``relationship: part_of`` edges, matching
    the propagation semantics of clusterProfiler::buildGOmap (which uses GO.db's
    GO{BP,MF,CC}ANCESTOR tables, built from is_a + part_of). ``regulates`` and
    its variants are intentionally excluded -- GO.db's ANCESTOR tables do not
    include them, and go-basic already drops cross-ontology edges.
    """
    terms = {}
    parents = defaultdict(set)
    in_term = False
    current_id = None
    current_name = None
    current_ns = None

    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line == "[Term]":
                if current_id and current_ns:
                    short_ns = NS_MAP.get(current_ns, current_ns)
                    terms[current_id] = {"name": current_name or "", "namespace": short_ns}
                in_term = True
                current_id = None
                current_name = None
                current_ns = None
            elif line.startswith("[") and line.endswith("]"):
                if current_id and current_ns:
                    short_ns = NS_MAP.get(current_ns, current_ns)
                    terms[current_id] = {"name": current_name or "", "namespace": short_ns}
                in_term = False
                current_id = None
            elif in_term:
                if line.startswith("id: "):
                    current_id = line[4:].strip()
                elif line.startswith("name: "):
                    current_name = line[6:].strip()
                elif line.startswith("namespace: "):
                    current_ns = line[11:].strip()
                elif line.startswith("is_a: "):
                    parent_id = line[6:].split("!")[0].strip()
                    if current_id:
                        parents[current_id].add(parent_id)
                elif line.startswith("relationship: part_of "):
                    # part_of propagates like is_a for annotation (a gene in a
                    # part is a gene in the whole). Other relationship types
                    # (regulates, has_part, ...) are NOT ancestors -- skip them.
                    parent_id = line[len("relationship: part_of "):].split("!")[0].strip()
                    if current_id and parent_id.startswith("GO:"):
                        parents[current_id].add(parent_id)
                elif line.startswith("is_obsolete: true"):
                    in_term = False
                    current_id = None

    if current_id and current_ns:
        short_ns = NS_MAP.get(current_ns, current_ns)
        terms[current_id] = {"name": current_name or "", "namespace": short_ns}

    log_msg("GO OBO: ", len(terms), " terms parsed")
    return terms, dict(parents)


def get_ancestors(go_id, parents, _cache=None):
    """Get all ancestor GO IDs (transitive is_a + part_of closure)."""
    if _cache is None:
        _cache = {}
    if go_id in _cache:
        return _cache[go_id]
    result = set()
    for parent in parents.get(go_id, []):
        result.add(parent)
        result |= get_ancestors(parent, parents, _cache)
    result.discard("all")
    _cache[go_id] = result
    return result


def build_ancestor_cache(parents):
    """Pre-compute ancestors for all GO IDs."""
    cache = {}
    for go_id in parents:
        get_ancestors(go_id, parents, cache)
    return cache
