"""KEGG REST API downloads and parsers."""

import hashlib
import re
from collections import defaultdict

from .utils import cached_download, log_msg


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


def kegg_get_batched(prefix, ids, cache_dir, refresh=False, *, retries=0,
                     rate_limit_s=0.0):
    """Fetch KEGG flat records for *ids* via batched ``get`` (<=10 per request).

    *prefix* is ``"cpd"`` or ``"rn"``. ids are sorted+de-duplicated so batching
    is deterministic, and each batch is cached under a filename bound to that
    batch's contents (a short hash) -- so an identical id set is a cache hit,
    while a changed set re-fetches instead of silently reusing a stale batch.
    ``retries``/``rate_limit_s`` are passed through to :func:`cached_download`.
    Returns the concatenated raw text of all records (KEGG separates records
    with ``///`` lines, which are preserved).
    """
    ids = sorted({i for i in ids if i})
    parts = []
    for i in range(0, len(ids), KEGG_GET_MAX):
        chunk = ids[i:i + KEGG_GET_MAX]
        url = "https://rest.kegg.jp/get/" + "+".join(f"{prefix}:{x}" for x in chunk)
        key = hashlib.sha1("|".join(chunk).encode()).hexdigest()[:12]
        path = cached_download(url, f"kegg_get_{prefix}_{key}.txt", cache_dir,
                               refresh, retries=retries, rate_limit_s=rate_limit_s)
        with open(path) as f:
            parts.append(f.read())
    return "".join(parts)


KEGG_REACTION_LIST_URL = "https://rest.kegg.jp/list/reaction"
_RXN_ID = re.compile(r"^R\d{5}$")


def download_reaction_list(cache_dir, refresh=False):
    """list/reaction: every KEGG reaction id (+ name), one cached file."""
    return cached_download(KEGG_REACTION_LIST_URL, "kegg_reaction_list.txt",
                           cache_dir, refresh)


def parse_reaction_list(path):
    """Parse list/reaction -> sorted unique reaction ids (``R#####``)."""
    ids = set()
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            rid = line.split("\t", 1)[0].replace("rn:", "").strip()
            if _RXN_ID.match(rid):
                ids.add(rid)
    return sorted(ids)


def fetch_and_reconcile_reactions(reaction_ids, cache_dir, *, refresh=False,
                                  retries=3, rate_limit_s=0.34):
    """Fetch + parse KEGG reaction records, reconciling requested vs parsed ids.

    Success is defined ONLY by reconciliation: a reaction id counts as retrieved
    when a record for it was actually parsed, never merely because a batch cache
    file exists. Unresolved ids are retried with the cache bypassed
    (``refresh=True``) -- the final pass fetches each missing id on its own, so a
    single bad/obsolete id in a batch cannot mask the rest. Returns
    ``(parsed, missing)``: *parsed* maps reaction id -> record, *missing* is the
    sorted list still unretrieved after all retries.
    """
    requested = sorted({r for r in reaction_ids if r})
    parsed = {}

    def absorb(text):
        for rec in iter_kegg_records(text):
            r = parse_reaction_record(rec)
            if r["id"] and r["id"] not in parsed:
                parsed[r["id"]] = r

    def missing_now():
        return [r for r in requested if r not in parsed]

    # Pass 0: normal batched fetch (cache allowed). Guarded so one bad id in a
    # chunk cannot abort the whole run -- the retry passes isolate it.
    if requested:
        try:
            absorb(kegg_get_batched("rn", requested, cache_dir, refresh=refresh,
                                    retries=1, rate_limit_s=rate_limit_s))
        except Exception as exc:  # noqa: BLE001 - reconciled via retries below
            log_msg("  initial reaction batch error (", exc,
                    "); reconciling via retries")

    missing = missing_now()
    for attempt in range(1, retries + 1):
        if not missing:
            break
        # Final pass isolates each id (smallest possible batch) with the cache
        # bypassed, so an incomplete/empty earlier response is never reused.
        per_id = attempt == retries
        if per_id:
            for rid in missing:
                try:
                    absorb(kegg_get_batched("rn", [rid], cache_dir, refresh=True,
                                            retries=1, rate_limit_s=rate_limit_s))
                except Exception:  # noqa: BLE001 - stays in `missing` -> fetch_failed
                    pass
        else:
            try:
                absorb(kegg_get_batched("rn", missing, cache_dir, refresh=True,
                                        retries=1, rate_limit_s=rate_limit_s))
            except Exception:  # noqa: BLE001 - isolated on the per-id pass
                pass
        missing = missing_now()

    return parsed, missing


def download_kegg_org_ko_links(kegg_code, cache_dir, refresh=False):
    """link/ko/<org>: organism gene -> KO pairs.

    KEGG does NOT link organism genes directly to reactions
    (``link/reaction/<org>`` is an invalid query, HTTP 400). Genes map to KOs,
    and KOs map to reactions KEGG-wide -- so the organism's reaction set is
    resolved through KO (see :func:`download_ko_reaction_links`).
    """
    url = f"https://rest.kegg.jp/link/ko/{kegg_code}"
    return cached_download(url, f"kegg_{kegg_code}_gene2ko.txt", cache_dir, refresh)


def parse_org_ko_links(path):
    """Parse link/ko/<org> -> (ordered KO ids, gene->KOs dict)."""
    ko_ids, gene2ko = [], defaultdict(list)
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            gene = re.sub(r"^[^:]+:", "", parts[0])
            ko = parts[1].replace("ko:", "")
            gene2ko[gene].append(ko)
            ko_ids.append(ko)
    return _dedup(ko_ids), dict(gene2ko)


def download_ko_reaction_links(cache_dir, refresh=False):
    """link/reaction/ko: KO -> reaction pairs, KEGG-wide (one cached file)."""
    url = "https://rest.kegg.jp/link/reaction/ko"
    return cached_download(url, "kegg_ko2reaction.txt", cache_dir, refresh)


def parse_ko_reaction_links(path):
    """Parse link/reaction/ko -> {KO: [reaction ids]} (KEGG-wide)."""
    ko2rxn = defaultdict(list)
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            ko = parts[0].replace("ko:", "")
            rxn = parts[1].replace("rn:", "")
            ko2rxn[ko].append(rxn)
    return {ko: _dedup(rxns) for ko, rxns in ko2rxn.items()}


def download_kegg_org_pathways(kegg_code, cache_dir, refresh=False):
    """list/pathway/<org>: organism pathway ids + names.

    Pass ``kegg_code=None`` for the KEGG-wide reference pathway list
    (``list/pathway`` -> ``map#####``), used by the KAAS (non-model) path.
    """
    suffix = f"/{kegg_code}" if kegg_code else ""
    dest = f"kegg_{kegg_code}_pathways.txt" if kegg_code else "kegg_ref_pathways.txt"
    return cached_download(f"https://rest.kegg.jp/list/pathway{suffix}", dest,
                           cache_dir, refresh)


def parse_org_pathways(path):
    """Parse list/pathway[/<org>] -> {pathway_id: name}, dropping any org suffix."""
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
    """info/<target>: used (best-effort) to record the KEGG snapshot dates."""
    url = f"https://rest.kegg.jp/info/{target}"
    return cached_download(url, f"kegg_info_{target}.txt", cache_dir, refresh)


# db line, e.g. "compound    19,584  2026/07/06"  (no release number anymore)
_KEGG_INFO_DB = re.compile(r"^\s*([a-z_]+)\s+[\d,]+\s+(\d{4})/(\d{2})/(\d{2})\s*$")


def parse_kegg_info_dates(path):
    """Parse info/kegg -> {db_name: 'YYYY-MM-DD'} from the per-database lines.

    Current info/kegg carries no release number; each database line is
    ``<name>   <count>   <YYYY/MM/DD>`` (pathway, ko, compound, reaction, ...).
    """
    dates = {}
    with open(path) as f:
        for line in f:
            m = _KEGG_INFO_DB.match(line)
            if m:
                dates[m.group(1)] = f"{m.group(2)}-{m.group(3)}-{m.group(4)}"
    return dates


def kegg_snapshot_version(dates, dbs=("pathway", "reaction", "compound")):
    """'KEGG snapshot <newest date>' across *dbs* (the databases this model is
    built from), or None if no dates were found."""
    picked = [dates[d] for d in dbs if d in dates] or list(dates.values())
    return f"KEGG snapshot {max(picked)}" if picked else None


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
    # The only reaction arrow KEGG uses in EQUATION is "<=>". Record the literal
    # token so a consumer never has to re-guess it, and leave it None when the
    # equation is absent or uses an unrecognized arrow (both are reported, not
    # silently dropped, downstream).
    arrow = "<=>" if "<=>" in equation else None
    reactants, products = [], []
    if arrow:
        left, right = equation.split(arrow, 1)
        reactants = _dedup(_CPD_ID.findall(left))
        products = _dedup(_CPD_ID.findall(right))
    enzymes = _dedup(_EC.findall(" ".join(rec.get("ENZYME", []))))
    ref_pathways = _dedup(_REF_PATH.findall(" ".join(rec.get("PATHWAY", []))))
    return {"id": rid, "name": name, "reactants": reactants, "products": products,
            "enzymes": enzymes, "ref_pathways": ref_pathways,
            # equation is the NORMALIZED full EQUATION field (KEGG wraps long
            # equations across lines; iter_kegg_records rejoins them with spaces)
            # -- not byte-verbatim source text.
            "equation": equation, "equation_arrow": arrow}
