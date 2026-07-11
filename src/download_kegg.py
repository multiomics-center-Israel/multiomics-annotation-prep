"""KEGG REST API downloads and parsers."""

import hashlib
import re
from collections import defaultdict

from .utils import cached_download


KEGG_KO_LIST_URL = "https://rest.kegg.jp/list/ko"
KEGG_KO_PATHWAY_URL = "https://rest.kegg.jp/link/pathway/ko"
KEGG_PATHWAY_LIST_URL = "https://rest.kegg.jp/list/pathway"


def download_kegg_rest(cache_dir, refresh=False):
    ko_name_path = cached_download(KEGG_KO_LIST_URL, "kegg_ko_to_name.txt",
                                   cache_dir, refresh)
    ko_path_path = cached_download(KEGG_KO_PATHWAY_URL, "kegg_ko_to_path.txt",
                                   cache_dir, refresh)
    pw_name_path = cached_download(KEGG_PATHWAY_LIST_URL, "kegg_pathway_names.txt",
                                   cache_dir, refresh)
    return ko_name_path, ko_path_path, pw_name_path


def parse_ko_to_name(path):
    """Parse ko -> {names, title, ec} from KEGG REST list/ko output."""
    result = {}
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t", 1)
            if len(parts) < 2:
                continue
            ko_id = parts[0].replace("ko:", "")
            desc = parts[1]
            names = desc
            title = ""
            ec = ""
            if "; " in desc:
                names, title = desc.split("; ", 1)
            ec_match = re.search(r"\[EC:[^\]]*\]", title)
            if ec_match:
                ec = ec_match.group(0)
                title = title.replace(ec, "").strip()
            result[ko_id] = {"names": names, "title": title, "ec": ec}
    return result


def parse_ko_to_path(path):
    """Parse ko -> list of reference pathway IDs (map*)."""
    ko2path = defaultdict(list)
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t", 1)
            if len(parts) < 2:
                continue
            ko = parts[0].replace("ko:", "")
            pth = parts[1].replace("path:", "")
            if pth.startswith("map"):
                ko2path[ko].append(pth)
    return dict(ko2path)


def parse_pathway_names(path):
    """Parse pathway ID -> pathway name."""
    result = {}
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t", 1)
            if len(parts) < 2:
                continue
            pw_id = parts[0].replace("path:", "")
            result[pw_id] = parts[1]
    return result


# ---------------------------------------------------------------------------
# Compound / reaction / pathway entities for the mummichog metabolic model.
#
# These are additive to the gene-set downloads above and reuse cached_download()
# so runs stay reproducible and offline after the first fetch. The mummichog
# model is COMPOUND-centric (compounds, real substrate/product reaction links,
# pathways-of-reactions) and is pulled directly from KEGG entities -- it is NOT
# derived from the gene sets above.
# ---------------------------------------------------------------------------

KEGG_GET_MAX = 10  # KEGG REST 'get' accepts up to 10 entries per request

_CPD_ID = re.compile(r"C\d{5}")
_EC = re.compile(r"\d+\.\d+\.\d+\.[0-9-]+")
_REF_PATH = re.compile(r"(?:rn|map)(\d{5})")


def _dedup(seq):
    """De-duplicate preserving first-seen order."""
    seen, out = set(), []
    for x in seq:
        if x not in seen:
            seen.add(x)
            out.append(x)
    return out


def kegg_get_batched(prefix, ids, cache_dir, refresh=False):
    """Fetch KEGG flat records for *ids* via batched ``get`` (<=10 per request).

    *prefix* is ``"cpd"`` or ``"rn"``. ids are sorted+de-duplicated so batching
    is deterministic, and each batch is cached under a filename bound to that
    batch's contents (a short hash) -- so an identical id set is a cache hit,
    while a changed set re-fetches instead of silently reusing a stale batch.
    Returns the concatenated raw text of all records (KEGG separates records
    with ``///`` lines, which are preserved).
    """
    ids = sorted({i for i in ids if i})
    parts = []
    for i in range(0, len(ids), KEGG_GET_MAX):
        chunk = ids[i:i + KEGG_GET_MAX]
        url = "https://rest.kegg.jp/get/" + "+".join(f"{prefix}:{x}" for x in chunk)
        key = hashlib.sha1("|".join(chunk).encode()).hexdigest()[:12]
        path = cached_download(url, f"kegg_get_{prefix}_{key}.txt", cache_dir, refresh)
        with open(path) as f:
            parts.append(f.read())
    return "".join(parts)


def download_kegg_org_reaction_links(kegg_code, cache_dir, refresh=False):
    """link/reaction/<org>: organism gene -> reaction pairs."""
    url = f"https://rest.kegg.jp/link/reaction/{kegg_code}"
    return cached_download(url, f"kegg_{kegg_code}_gene2reaction.txt",
                           cache_dir, refresh)


def parse_org_reaction_links(path):
    """Parse link/reaction/<org> -> (ordered reaction ids, gene->reactions dict).

    The reaction set defines which reactions are present in the organism, keyed
    on the organism's own genes (real enzyme presence).
    """
    reaction_ids, gene2rxn = [], defaultdict(list)
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            gene = re.sub(r"^[^:]+:", "", parts[0])
            rxn = parts[1].replace("rn:", "")
            gene2rxn[gene].append(rxn)
            reaction_ids.append(rxn)
    return _dedup(reaction_ids), dict(gene2rxn)


def download_kegg_org_pathways(kegg_code, cache_dir, refresh=False):
    """list/pathway/<org>: organism pathway ids + names."""
    url = f"https://rest.kegg.jp/list/pathway/{kegg_code}"
    return cached_download(url, f"kegg_{kegg_code}_pathways.txt", cache_dir, refresh)


def parse_org_pathways(path):
    """Parse list/pathway/<org> -> {pathway_id: name}, dropping the org suffix."""
    result = {}
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t", 1)
            if len(parts) < 2:
                continue
            pid = parts[0].replace("path:", "")
            name = parts[1]
            if " - " in name:  # "<name> - <Organism ...>"
                name = name.rsplit(" - ", 1)[0]
            result[pid] = name.strip()
    return result


def download_kegg_info(target, cache_dir, refresh=False):
    """info/<target>: used (best-effort) to record the KEGG release version."""
    url = f"https://rest.kegg.jp/info/{target}"
    return cached_download(url, f"kegg_info_{target}.txt", cache_dir, refresh)


def parse_kegg_release(path):
    """Extract a 'KEGG Release ...' string from an info/<target> file, or None."""
    with open(path) as f:
        for line in f:
            m = re.search(r"Release\s+(\d+\.\d+)", line)
            if m:
                return f"KEGG Release {m.group(1)}"
    return None


def iter_kegg_records(text):
    """Yield one dict per KEGG flat record: {FIELD: [value_line, ...]}.

    KEGG format puts the field label in columns 1-12 and the value from column
    13 on; continuation lines are blank in columns 1-12; records end with a
    line beginning ``///``.
    """
    record, current = {}, None
    for raw in text.splitlines():
        if raw.startswith("///"):
            if record:
                yield record
            record, current = {}, None
            continue
        if not raw.strip():
            continue
        label, value = raw[:12].strip(), raw[12:].rstrip()
        if label:
            current = label
            record.setdefault(current, [])
            if value:
                record[current].append(value)
        elif current is not None and value:
            record[current].append(value)
    if record:
        yield record


def parse_compound_record(rec):
    """Parse one KEGG compound record -> {id, name, formula, exact_mass}."""
    entry = rec.get("ENTRY", [""])[0].split()
    cid = entry[0] if entry else ""
    name = ""
    if rec.get("NAME"):
        name = rec["NAME"][0].strip().rstrip(";").strip()
    formula = rec.get("FORMULA", [""])[0].strip() if rec.get("FORMULA") else ""
    exact_mass = None
    if rec.get("EXACT_MASS"):
        try:
            exact_mass = float(rec["EXACT_MASS"][0].split()[0])
        except (ValueError, IndexError):
            exact_mass = None
    return {"id": cid, "name": name, "formula": formula, "exact_mass": exact_mass}


def parse_reaction_record(rec):
    """Parse one KEGG reaction record.

    Returns {id, name, reactants, products, enzymes, ref_pathways}. Reactants and
    products come from the EQUATION (left/right of ``<=>``) so the links are the
    real substrate/product directions, not pathway co-occurrence. ref_pathways
    are the 5-digit reference pathway numbers this reaction belongs to.
    """
    entry = rec.get("ENTRY", [""])[0].split()
    rid = entry[0] if entry else ""
    name = ""
    if rec.get("NAME"):
        name = rec["NAME"][0].strip().rstrip(";").strip()
    equation = " ".join(rec.get("EQUATION", []))
    reactants, products = [], []
    if "<=>" in equation:
        left, right = equation.split("<=>", 1)
        reactants = _dedup(_CPD_ID.findall(left))
        products = _dedup(_CPD_ID.findall(right))
    enzymes = _dedup(_EC.findall(" ".join(rec.get("ENZYME", []))))
    ref_pathways = _dedup(_REF_PATH.findall(" ".join(rec.get("PATHWAY", []))))
    return {"id": rid, "name": name, "reactants": reactants, "products": products,
            "enzymes": enzymes, "ref_pathways": ref_pathways}
