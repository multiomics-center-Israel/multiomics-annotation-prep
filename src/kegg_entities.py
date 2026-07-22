"""Organism KEGG source loading, shared by the compound-centric modules.

Both the mummichog metabolic model (:mod:`prepare_mummichog_model`) and the KEGG
compound-set GMT (:mod:`prepare_kegg_compound_sets`) are built from the SAME
source shape: a KO list -> reactions -> compounds -> pathways. That loading step
lives here so both modules reuse it without importing each other (which would be
circular). It is deliberately dependency-light: only stdlib + the
``requests``-backed downloads in :mod:`download_kegg`. The heavy, model-specific
deps (``metDataModel`` for the structure, ``mass2chem`` for the masses) live in
the modules that need them, imported lazily there.
"""

from collections import defaultdict

from .download_kegg import (
    _dedup,
    download_kegg_info,
    download_kegg_org_ko_links,
    download_kegg_org_pathways,
    download_ko_reaction_links,
    iter_kegg_records,
    kegg_get_batched,
    kegg_snapshot_version,
    parse_compound_record,
    parse_kegg_info_dates,
    parse_ko_reaction_links,
    parse_org_ko_links,
    parse_org_pathways,
    parse_reaction_record,
)
from .utils import log_msg

SOURCE_TOKENS = {"kegg_org": "kegg", "kaas": "kaas"}


def is_metabolic_map(number):
    """True for a real metabolic pathway map (KEGG ``00xxx``, i.e. 00001-00999).

    Excludes KEGG's "Global and overview maps" (the 01100-01299 band: Metabolic
    pathways, Biosynthesis of secondary metabolites, Carbon metabolism,
    Biosynthesis of amino acids, ...) and any non-metabolic category. Those
    overview maps aggregate many pathways and distort enrichment, so they are
    kept out of the model's / gene-set's pathway list.
    """
    try:
        return 0 < int(number) < 1000
    except (TypeError, ValueError):
        return False


# ---------------------------------------------------------------------------
# Source loading. Both the organism-code path and the KAAS path reduce to the
# SAME shape: a KO list -> reactions -> compounds -> pathways. The only thing
# that differs is where the KO list comes from.
#
# KEGG does not link organism genes directly to reactions
# (link/reaction/<org> is invalid). The organism's reactions are resolved
# through KO: link/ko/<org> (gene->KO) intersected with link/reaction/ko
# (KO->reaction, KEGG-wide).
# ---------------------------------------------------------------------------

def parse_kaas_ko_list(kaas_file):
    """Parse a KAAS query.ko.txt (gene<TAB>KO) -> (ordered KO ids, gene->KOs)."""
    ko_ids, gene2ko = [], defaultdict(list)
    with open(kaas_file) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t", 1)
            gene = parts[0].strip()
            ko = parts[1].strip().replace("ko:", "") if len(parts) > 1 else ""
            if not ko:
                continue
            gene2ko[gene].append(ko)
            ko_ids.append(ko)
    return _dedup(ko_ids), dict(gene2ko)


def resolve_ko_list(source, kegg_code, cache_dir, refresh, kaas_file):
    """The unifying step: get the organism's KO list from either source.

    Returns ``(ko_list, {"n_genes": ...})``.
    """
    if source == "kegg_org":
        if not kegg_code:
            raise ValueError("kegg_code is required for source='kegg_org'")
        log_msg("KO list from KEGG organism code: ", kegg_code)
        ko_list, gene2ko = parse_org_ko_links(
            download_kegg_org_ko_links(kegg_code, cache_dir, refresh))
    elif source == "kaas":
        if not kaas_file:
            raise ValueError("kaas_file is required for source='kaas'")
        log_msg("KO list from KAAS file: ", kaas_file)
        ko_list, gene2ko = parse_kaas_ko_list(kaas_file)
    else:
        raise ValueError(
            f"Unknown source {source!r} (expected 'kegg_org' or 'kaas')")
    log_msg("  genes: ", len(gene2ko), "; distinct KOs: ", len(ko_list))
    return ko_list, {"n_genes": len(gene2ko)}


def load_source(source, kegg_code, cache_dir, refresh=False, kaas_file=None):
    """KO list -> reactions -> compounds -> pathways (source-agnostic core).

    Returns ``{reactions, compounds, pathway_names, pathway_prefix, source,
    source_version, ko_coverage}``.
    """
    ko_list, ko_meta = resolve_ko_list(source, kegg_code, cache_dir, refresh,
                                       kaas_file)

    # KO -> reaction (KEGG-wide), intersected with the organism's KOs.
    ko2rxn = parse_ko_reaction_links(download_ko_reaction_links(cache_dir, refresh))
    kos_with_rxn = [ko for ko in ko_list if ko2rxn.get(ko)]
    reaction_ids = _dedup([r for ko in ko_list for r in ko2rxn.get(ko, [])])
    log_msg("  KOs mapping to >=1 reaction: ", len(kos_with_rxn), "/",
            len(ko_list), "; reactions: ", len(reaction_ids))

    rn_text = kegg_get_batched("rn", reaction_ids, cache_dir, refresh)
    reactions = [r for r in (parse_reaction_record(rec)
                             for rec in iter_kegg_records(rn_text)) if r["id"]]

    cpd_ids = _dedup([c for r in reactions
                      for c in (r["reactants"] + r["products"])])
    log_msg("  distinct compounds referenced by reactions: ", len(cpd_ids))
    cpd_text = kegg_get_batched("cpd", cpd_ids, cache_dir, refresh)
    compounds = {}
    for rec in iter_kegg_records(cpd_text):
        c = parse_compound_record(rec)
        if c["id"]:
            compounds[c["id"]] = c

    # Pathways: organism pathways for a KEGG code (cre#####); reference maps
    # (map#####) for KAAS, since a non-model organism has none of its own.
    if source == "kegg_org":
        pathway_names = parse_org_pathways(
            download_kegg_org_pathways(kegg_code, cache_dir, refresh))
        pathway_prefix, source_str = kegg_code, "KEGG REST"
    else:
        pathway_names = parse_org_pathways(
            download_kegg_org_pathways(None, cache_dir, refresh))
        pathway_prefix, source_str = "map", "KEGG REST (KAAS KO list)"
    log_msg("  candidate pathways: ", len(pathway_names))

    # info/kegg no longer has a release number, only per-database snapshot dates;
    # provenance is the newest date across the dbs this model is built from.
    source_version, kegg_db_dates = None, {}
    try:
        all_dates = parse_kegg_info_dates(
            download_kegg_info("kegg", cache_dir, refresh))
        source_version = kegg_snapshot_version(all_dates)
        kegg_db_dates = {db: all_dates.get(db)
                         for db in ("pathway", "reaction", "compound")}
    except Exception as exc:  # noqa: BLE001 - version is best-effort, not fatal
        log_msg("  (KEGG snapshot date unavailable: ", exc, ")")

    return {
        "reactions": reactions,
        "compounds": compounds,
        "pathway_names": pathway_names,
        "pathway_prefix": pathway_prefix,
        "source": source_str,
        "source_version": source_version,
        "kegg_db_dates": kegg_db_dates,
        "ko_coverage": {
            "n_genes": ko_meta.get("n_genes"),
            "n_kos": len(ko_list),
            "n_kos_with_reaction": len(kos_with_rxn),
            "n_reactions_from_kos": len(reaction_ids),
        },
    }
