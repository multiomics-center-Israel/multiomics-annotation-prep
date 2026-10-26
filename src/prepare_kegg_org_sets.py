"""KEGG pathway gene sets for a KEGG organism code, re-keyed to project ids.

For organisms KEGG already annotates (e.g. Giardia, ``gla``) there is no need
to run KAAS: KEGG's own gene -> pathway links are fetched over REST. KEGG keys
them on its organism gene ids (for many genomes the NCBI locus tags), so an
optional id map (``id_maps``) re-keys them onto the ids in the project's data
-- the same map the GO side uses, so GO and KEGG sets share one id space.
Writes the KEGG .tab tables, ``KEGG_pathway.gmt`` and ``KEGG_provenance.json``.
"""

import json
import os
import re
from collections import defaultdict
from datetime import datetime, timezone

from .download_kegg_org import kegg_org_links
from .id_maps import resolve_id_map
from .utils import (cached_download, git_head_sha, log_msg, sha256_file, write_gmt,
                    write_pathway2gene, write_pathway2name)


def kegg_release_info(kegg_org, cache_dir, refresh=False):
    """Return the non-empty lines of ``/info/<org>`` (database name + release)."""
    path = cached_download(f"https://rest.kegg.jp/info/{kegg_org}",
                           f"kegg_{kegg_org}_info.txt", cache_dir, refresh)
    with open(path) as f:
        return [line.strip() for line in f if line.strip()][:10]


def remap_pathway_genes(pathways, genes, mapping):
    """Re-key ``(pathway, gene)`` links onto the targets of *mapping*.

    A gene with several targets gives the pathway each of them; unmapped genes
    are dropped. Raises when nothing maps, since that means the map's source
    ids are not the KEGG gene ids.
    """
    pairs = []
    unmapped = set()
    for p, g in zip(pathways, genes):
        targets = mapping.get(g)
        if not targets:
            unmapped.add(g)
            continue
        pairs.extend((p, t) for t in targets)
    if not pairs:
        raise RuntimeError("id map matched none of the KEGG gene ids -- check that "
                           "its source column uses the KEGG organism's gene ids")
    if unmapped:
        log_msg("id map: ", len(unmapped), " of ", len(set(genes)),
                " KEGG genes have no mapping and were dropped")
    return pairs


def prepare_kegg_org_sets(kegg_org, out_dir, cache_dir, id_map=None,
                          ncbi_assembly=None, refresh=False):
    """Write KEGG tables, ``KEGG_pathway.gmt`` and provenance for an organism.

    kegg_org:      KEGG organism code (e.g. ``gla``).
    id_map / ncbi_assembly: optional re-keying, as in ``id_maps.resolve_id_map``.
    """
    if not re.match(r"^[a-z]{3,4}$", kegg_org):
        raise ValueError(f"'{kegg_org}' is not a KEGG organism code (e.g. gla, hsa)")
    genes, pathways, pw_names = kegg_org_links(kegg_org, cache_dir, refresh)
    if not genes:
        raise RuntimeError(f"KEGG returned no gene -> pathway links for '{kegg_org}'")
    genes = [re.sub(r"^[^:]+:", "", g) for g in genes]
    log_msg("KEGG (org=", kegg_org, "): ", len(set(pathways)), " pathways, ",
            len(set(genes)), " genes")

    mapping, map_info = resolve_id_map(id_map, ncbi_assembly, cache_dir, refresh)
    if mapping is not None:
        pairs = remap_pathway_genes(pathways, genes, mapping)
    else:
        pairs = list(zip(pathways, genes))

    write_pathway2gene(pairs, os.path.join(out_dir, "KEGG_pathway2gene.tab"))
    used = sorted({p for p, _ in pairs})
    write_pathway2name([(p, pw_names.get(p, "")) for p in used],
                       os.path.join(out_dir, "KEGG_pathway2name.tab"))
    term2genes = defaultdict(list)
    for p, g in pairs:
        term2genes[p].append(g)
    write_gmt(dict(term2genes), pw_names, os.path.join(out_dir, "KEGG_pathway.gmt"))

    links = os.path.join(cache_dir, f"kegg_{kegg_org}_gene2path.txt")
    provenance = {
        "kegg_source": "KEGG REST",
        "kegg_org": kegg_org,
        "kegg_release": kegg_release_info(kegg_org, cache_dir, refresh),
        "links": {"url": f"https://rest.kegg.jp/link/pathway/{kegg_org}",
                  "sha256": sha256_file(links)},
        "id_map": map_info,
        "pathways": len(used),
        "annotated_ids": len({g for _, g in pairs}),
        "builder_git_sha": git_head_sha() or "unknown",
        "build_timestamp_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
    }
    path = os.path.join(out_dir, "KEGG_provenance.json")
    with open(path, "w") as f:
        json.dump(provenance, f, indent=2)
    log_msg("wrote provenance: ", path)
