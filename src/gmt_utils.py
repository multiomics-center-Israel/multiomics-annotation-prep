"""GMT gene-set helpers — Python port of R/gmt_utils.R.

Read/write/validate/filter/merge for GMT files, plus KEGG- and biomaRt-backed
gene-set generation. Mirrors the vendored R helpers so the Python and R paths
produce the same GMT files. Every remote fetch goes through cached_download()
so runs are reproducible and offline after the first fetch.

A gene set collection is a `GeneSets` (an ordered dict term_id -> [members])
carrying a parallel `.descriptions` dict term_id -> name. This mirrors the R
convention of a named list with a "descriptions" attribute.
"""

import hashlib
import os
import re
import statistics

from .utils import (cached_download, cached_download_filtered, cached_post,
                    ensure_dir, log_msg)

KEGG_REST = "https://rest.kegg.jp"
BIOMART_HOST = "https://www.ensembl.org"
REACTOME_DOWNLOAD = "https://reactome.org/download/current"


class GeneSets(dict):
    """Ordered term_id -> [members] mapping with a parallel descriptions dict."""

    def __init__(self, *args, descriptions=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.descriptions = dict(descriptions or {})

    def _subset_descriptions(self):
        """Restrict descriptions to the currently present term ids."""
        self.descriptions = {k: self.descriptions.get(k, "") for k in self}
        return self


# ==============================================================================
# GMT FILE READING / WRITING
# ==============================================================================

def read_gmt(gmt_file):
    """Load a GMT file into a GeneSets (term_id -> members, plus descriptions)."""
    if not os.path.exists(gmt_file):
        raise FileNotFoundError(f"GMT file not found: {gmt_file}")

    sets = GeneSets()
    with open(gmt_file) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                name = parts[0]
                sets[name] = [g for g in parts[2:] if g != ""]
                sets.descriptions[name] = parts[1]

    log_msg(f"Loaded {len(sets)} gene sets from GMT file")
    return sets


def write_gmt(pathways, output_file, descriptions=None, verbose=True):
    """Write a GeneSets to a GMT file (term_id<tab>description<tab>members...)."""
    if descriptions is None:
        descriptions = getattr(pathways, "descriptions", {})

    ensure_dir(os.path.dirname(output_file) or ".")
    with open(output_file, "w") as f:
        for name, genes in pathways.items():
            desc = descriptions.get(name, "") or ""
            f.write("\t".join([name, desc, *genes]) + "\n")

    if verbose:
        log_msg(f"Wrote {len(pathways)} pathways to: {output_file}")
    return output_file


# ==============================================================================
# GMT VALIDATION / MANIPULATION
# ==============================================================================

def validate_gmt(gmt, feature_ids, min_coverage=0.1, min_pathway_size=5,
                 max_pathway_size=500, verbose=True):
    """Validate a GMT against a dataset's feature IDs; report coverage + sizes."""
    pathways = read_gmt(gmt) if isinstance(gmt, str) else gmt
    feature_set = set(feature_ids)

    all_genes = set()
    for genes in pathways.values():
        all_genes.update(genes)

    matching = all_genes & feature_set
    coverage = len(matching) / len(feature_set) if feature_set else 0.0
    gmt_coverage = len(matching) / len(all_genes) if all_genes else 0.0

    filtered = GeneSets(descriptions=getattr(pathways, "descriptions", {}))
    for name, genes in pathways.items():
        kept = [g for g in genes if g in feature_set]
        if min_pathway_size <= len(kept) <= max_pathway_size:
            filtered[name] = kept
    filtered._subset_descriptions()

    result = {
        "valid": coverage >= min_coverage,
        "coverage": coverage,
        "gmt_coverage": gmt_coverage,
        "n_pathways_original": len(pathways),
        "n_pathways_filtered": len(filtered),
        "n_genes_in_gmt": len(all_genes),
        "n_genes_matching": len(matching),
        "n_features_in_data": len(feature_set),
        "filtered_pathways": filtered,
        "matching_genes": matching,
    }

    if verbose:
        log_msg("=== GMT Validation Report ===")
        log_msg(f"GMT file: {result['n_pathways_original']} pathways, "
                f"{result['n_genes_in_gmt']} unique genes")
        log_msg(f"Your data: {result['n_features_in_data']} features")
        log_msg("Overlap:")
        log_msg(f"  {100 * coverage:.1f}% of your features found in GMT "
                f"({result['n_genes_matching']} genes)")
        log_msg(f"  {100 * gmt_coverage:.1f}% of GMT genes found in your data")
        log_msg(f"After filtering (size {min_pathway_size}-{max_pathway_size}):")
        log_msg(f"  {result['n_pathways_filtered']} pathways retained "
                f"(from {result['n_pathways_original']} original)")
        if not result["valid"]:
            log_msg(f"WARNING: Coverage ({100 * coverage:.1f}%) is below "
                    f"threshold ({100 * min_coverage:.1f}%)")
        else:
            log_msg(f"Validation PASSED: Coverage ({100 * coverage:.1f}%) "
                    f"meets threshold")

    return result


def filter_gmt_by_size(pathways, min_size=10, max_size=500):
    """Keep only sets whose member count is within [min_size, max_size]."""
    out = GeneSets(descriptions=getattr(pathways, "descriptions", {}))
    for name, genes in pathways.items():
        if min_size <= len(genes) <= max_size:
            out[name] = genes
    return out._subset_descriptions()


def merge_gmt(gmt_files, prefixes=None, verbose=True):
    """Merge several GMTs (paths or GeneSets), optionally prefixing set names."""
    merged = GeneSets()
    for i, src in enumerate(gmt_files):
        pathways = read_gmt(src) if isinstance(src, str) else src
        prefix = prefixes[i] if prefixes and i < len(prefixes) else None
        for name, genes in pathways.items():
            key = f"{prefix}_{name}" if prefix else name
            merged[key] = genes
            merged.descriptions[key] = pathways.descriptions.get(name, "") \
                if isinstance(pathways, GeneSets) else ""

    if verbose:
        log_msg(f"Merged {len(gmt_files)} GMT sources into {len(merged)} pathways")
    return merged


def convert_gmt_ids(pathways, mapping, verbose=True):
    """Remap member IDs through `mapping` (dict from_id -> to_id); drop unmapped."""
    converted = GeneSets(descriptions=getattr(pathways, "descriptions", {}))
    for name, genes in pathways.items():
        new_genes = list(dict.fromkeys(
            mapping[g] for g in genes if g in mapping and mapping[g] is not None))
        if new_genes:
            converted[name] = new_genes
    converted._subset_descriptions()

    if verbose:
        log_msg(f"Converted GMT: {len(pathways)} -> {len(converted)} pathways "
                f"(removed {len(pathways) - len(converted)} empty)")
    return converted


def gmt_summary(pathways):
    """Return a list of {pathway_id, description, n_genes} dicts."""
    descriptions = getattr(pathways, "descriptions", {})
    return [
        {"pathway_id": name,
         "description": descriptions.get(name, ""),
         "n_genes": len(genes)}
        for name, genes in pathways.items()
    ]


def _log_sizes(pathways):
    sizes = [len(g) for g in pathways.values()]
    if sizes:
        log_msg(f"Pathway sizes: min={min(sizes)}, "
                f"median={int(statistics.median(sizes))}, max={max(sizes)}")


# ==============================================================================
# KEGG REST HELPERS (replace KEGGREST)
# ==============================================================================

def kegg_list_pathways(organism, cache_dir, refresh=False):
    """Return dict pathway_id -> name for a KEGG organism code."""
    url = f"{KEGG_REST}/list/pathway/{organism}"
    path = cached_download(url, f"kegg_pathways_{organism}.txt", cache_dir, refresh)
    result = {}
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            pid, _, name = line.partition("\t")
            result[pid] = name
    return result


def kegg_get(ids, cache_dir, refresh=False):
    """Fetch KEGG flat-file entries for up to 10 ids; return parsed section dicts."""
    query = "+".join(ids)
    digest = hashlib.md5(query.encode()).hexdigest()[:12]
    path = cached_download(f"{KEGG_REST}/get/{query}",
                           f"kegg_get_{digest}.txt", cache_dir, refresh)
    with open(path) as f:
        return _parse_kegg_entries(f.read())


def _parse_kegg_entries(content):
    """Split a KEGG flat file into per-entry {SECTION: [lines]} dicts."""
    entries = []
    for block in content.split("\n///"):
        if not block.strip():
            continue
        sections = {}
        current = None
        for line in block.split("\n"):
            if not line.strip():
                continue
            key = line[:12].strip()
            body = line[12:].rstrip()
            if key:
                current = key
                sections.setdefault(current, []).append(body)
            elif current:
                sections[current].append(body)
        entries.append(sections)
    return entries


def parse_kegg_id_block(lines):
    """Parse a GENE/COMPOUND section into an ordered dict id -> description."""
    out = {}
    for line in lines:
        line = line.strip()
        if not line:
            continue
        parts = line.split(None, 1)
        out[parts[0]] = parts[1] if len(parts) > 1 else ""
    return out


# ==============================================================================
# GMT GENERATION FROM DATABASES
# ==============================================================================

def generate_gmt_from_kegg(organism, id_type="entrez", output_file=None,
                           cache_dir="data", refresh=False,
                           name_strip=r" - .*$", verbose=True):
    """Build a KEGG pathway -> genes GeneSets for an organism code.

    id_type "entrez" keeps the KEGG gene IDs; "symbol" takes the first token of
    each gene's description. `name_strip` is applied to each pathway name.
    """
    pathway_list = kegg_list_pathways(organism, cache_dir, refresh)
    if not pathway_list:
        raise ValueError(f"No KEGG pathways found for organism: {organism}")
    if verbose:
        log_msg(f"Found {len(pathway_list)} pathways, retrieving genes...")

    pathways = GeneSets()
    for i, (pid, pname) in enumerate(pathway_list.items(), start=1):
        pname = re.sub(name_strip, "", pname)
        try:
            entry = kegg_get([pid], cache_dir, refresh)[0]
        except Exception:
            if verbose:
                log_msg(f"  Warning: Failed to fetch {pid}")
            continue

        gene_block = parse_kegg_id_block(entry.get("GENE", []))
        if not gene_block:
            continue
        if id_type == "symbol":
            genes = [desc.split(";")[0].strip() for desc in gene_block.values()]
        else:
            genes = list(gene_block.keys())

        clean_id = pid.replace("path:", "")
        pathways[clean_id] = list(dict.fromkeys(g for g in genes if g))
        pathways.descriptions[clean_id] = pname

        if verbose and i % 25 == 0:
            log_msg(f"  Progress: {i}/{len(pathway_list)} pathways")

    if verbose:
        log_msg(f"Generated {len(pathways)} KEGG pathways")
        _log_sizes(pathways)
    if output_file:
        write_gmt(pathways, output_file, verbose=verbose)
    return pathways


# Ensembl dataset lookup, mirrored from R/gmt_utils.R (.gmt_ensembl_dataset).
_ENSEMBL_DATASETS = {
    "homo_sapiens": "hsapiens_gene_ensembl", "human": "hsapiens_gene_ensembl",
    "mus_musculus": "mmusculus_gene_ensembl", "mouse": "mmusculus_gene_ensembl",
    "rattus_norvegicus": "rnorvegicus_gene_ensembl", "rat": "rnorvegicus_gene_ensembl",
    "danio_rerio": "drerio_gene_ensembl", "zebrafish": "drerio_gene_ensembl",
    "drosophila_melanogaster": "dmelanogaster_gene_ensembl", "fly": "dmelanogaster_gene_ensembl",
    "caenorhabditis_elegans": "celegans_gene_ensembl", "worm": "celegans_gene_ensembl",
    "saccharomyces_cerevisiae": "scerevisiae_gene_ensembl", "yeast": "scerevisiae_gene_ensembl",
    "arabidopsis_thaliana": "athaliana_gene_ensembl",
    "gallus_gallus": "ggallus_gene_ensembl", "chicken": "ggallus_gene_ensembl",
    "sus_scrofa": "sscrofa_gene_ensembl", "pig": "sscrofa_gene_ensembl",
    "bos_taurus": "btaurus_gene_ensembl", "cow": "btaurus_gene_ensembl",
}

_GO_NAMESPACE = {
    "BP": "biological_process",
    "MF": "molecular_function",
    "CC": "cellular_component",
}


def ensembl_dataset(organism):
    """Return the Ensembl mart dataset for an organism name, or None."""
    key = re.sub(r"\s+", "_", organism.strip().lower())
    dataset = _ENSEMBL_DATASETS.get(key)
    if dataset is None:
        parts = key.split("_")
        if len(parts) >= 2:
            dataset = f"{parts[0][0]}{parts[1]}_gene_ensembl"
    return dataset


def biomart_query(dataset_name, attributes, cache_dir="data", refresh=False,
                  host=BIOMART_HOST):
    """Run a BioMart attribute query; return rows as lists of strings.

    Talks to the martservice REST endpoint directly rather than through
    `pybiomart`, whose per-dataset configuration XML fetch is large (~2 MB) and
    fails to parse whenever Ensembl truncates the response.
    """
    attrs = "".join(f'<Attribute name="{a}"/>' for a in attributes)
    query = (
        '<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query>'
        '<Query virtualSchemaName="default" formatter="TSV" header="0" '
        'uniqueRows="1" count="" datasetConfigVersion="0.6">'
        f'<Dataset name="{dataset_name}" interface="default">{attrs}</Dataset>'
        "</Query>"
    )
    digest = hashlib.md5(query.encode()).hexdigest()[:12]
    path = cached_post(f"{host}/biomart/martservice", {"query": query},
                       f"biomart_{dataset_name}_{digest}.tsv", cache_dir, refresh)

    with open(path, encoding="utf-8", errors="replace") as f:
        head = f.read(400)
        if "Query ERROR" in head or head.lstrip().startswith("<"):
            os.remove(path)  # don't cache an error page
            raise RuntimeError(
                f"BioMart returned an error for dataset '{dataset_name}':\n"
                f"{head.strip()[:300]}")
        f.seek(0)
        rows = [line.rstrip("\n").split("\t") for line in f if line.strip()]

    ncol = len(attributes)
    return [r for r in rows if len(r) == ncol]


def generate_gmt_from_go(organism, id_type="symbol", ontology="BP",
                         output_file=None, min_genes=10, max_genes=500,
                         cache_dir="data", refresh=False, verbose=True):
    """Build a GO -> genes GeneSets from Ensembl BioMart.

    Ensembl only covers model organisms; for non-model organisms use the
    GAF/eggNOG builders (build_ankri_go_gmt.py / build_elad_go_gmt.py) instead.
    """
    dataset_name = ensembl_dataset(organism)
    if dataset_name is None:
        raise ValueError(f"Could not determine Ensembl dataset for: {organism}")

    gene_attr = {
        "ensembl": "ensembl_gene_id",
        "symbol": "external_gene_name",
        "entrez": "entrezgene_id",
    }.get(id_type, "external_gene_name")

    if verbose:
        log_msg(f"Querying Ensembl BioMart (dataset: {dataset_name})...")
    rows = biomart_query(dataset_name,
                         [gene_attr, "go_id", "name_1006", "namespace_1003"],
                         cache_dir=cache_dir, refresh=refresh)
    if verbose:
        log_msg(f"Retrieved {len(rows)} gene-to-GO rows")

    wanted_ns = None if ontology == "all" else _GO_NAMESPACE.get(ontology, ontology)

    pathways = GeneSets()
    for gene, go_id, name, namespace in rows:
        if not gene or not go_id:
            continue
        if wanted_ns is not None and namespace != wanted_ns:
            continue
        pathways.setdefault(go_id, [])
        pathways[go_id].append(gene)
        pathways.descriptions.setdefault(go_id, name)
    for go_id, genes in pathways.items():
        pathways[go_id] = list(dict.fromkeys(genes))

    if verbose:
        log_msg(f"Generated {len(pathways)} GO pathways")
        _log_sizes(pathways)
    if min_genes is not None:
        pathways = filter_gmt_by_size(pathways, min_genes, max_genes)
    if output_file:
        write_gmt(pathways, output_file, verbose=verbose)
    return pathways


# ==============================================================================
# REACTOME
# ==============================================================================

def _is_gene_level(identifier):
    """True unless the ID is an Ensembl transcript (…T#) or protein (…P#) ID.

    Reactome's Ensembl mapping lists gene, transcript and protein IDs for the
    same annotation; only the gene rows should become GMT members. IDs that are
    not Ensembl-style are kept as-is.
    """
    m = re.match(r"^ENS[A-Z]*([GTP])\d+", identifier)
    return m is None or m.group(1) == "G"


def ensembl_symbol_map(organism, cache_dir="data", refresh=False):
    """Return dict ensembl_gene_id -> gene symbol for an organism, via BioMart."""
    dataset_name = ensembl_dataset(organism)
    if dataset_name is None:
        raise ValueError(f"Could not determine Ensembl dataset for: {organism}")
    rows = biomart_query(dataset_name, ["ensembl_gene_id", "external_gene_name"],
                         cache_dir=cache_dir, refresh=refresh)
    return {gid: sym for gid, sym in rows if gid and sym}


def generate_gmt_from_reactome(organism, id_type="symbol", output_file=None,
                               min_genes=10, max_genes=500, cache_dir="data",
                               refresh=False, verbose=True):
    """Build a Reactome pathway -> genes GeneSets for an organism.

    Uses Reactome's Ensembl2Reactome_All_Levels mapping (all levels of the
    pathway hierarchy). The download covers every species, so it is streamed and
    only the requested species' lines are cached. id_type "ensembl" keeps the
    Ensembl gene IDs; "symbol" maps them through BioMart.
    """
    species = organism.strip()
    species_tag = re.sub(r"\s+", "_", species.lower())
    suffix = f"\t{species}"

    path = cached_download_filtered(
        f"{REACTOME_DOWNLOAD}/Ensembl2Reactome_All_Levels.txt",
        f"reactome_ensembl_{species_tag}.txt", cache_dir,
        keep=lambda line: line.endswith(suffix), refresh=refresh)

    mapping = None
    if id_type == "symbol":
        if verbose:
            log_msg("Mapping Ensembl gene IDs to symbols...")
        mapping = ensembl_symbol_map(species, cache_dir, refresh)

    pathways = GeneSets()
    n_rows = n_dropped = 0
    with open(path, encoding="utf-8") as f:
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 4:
                continue
            gene, pathway_id, _, pathway_name = cols[:4]
            if not _is_gene_level(gene):
                continue
            n_rows += 1
            if mapping is not None:
                gene = mapping.get(gene)
                if gene is None:
                    n_dropped += 1
                    continue
            pathways.setdefault(pathway_id, [])
            pathways[pathway_id].append(gene)
            pathways.descriptions.setdefault(pathway_id, pathway_name)
    for pathway_id, genes in pathways.items():
        pathways[pathway_id] = list(dict.fromkeys(genes))

    if verbose:
        log_msg(f"Read {n_rows} gene-level Reactome rows for {species}")
        if n_dropped:
            log_msg(f"  {n_dropped} rows dropped: the Ensembl gene ID is retired "
                    f"or has no symbol in the current Ensembl release")
        log_msg(f"Generated {len(pathways)} Reactome pathways")
        _log_sizes(pathways)
    if min_genes is not None:
        pathways = filter_gmt_by_size(pathways, min_genes, max_genes)
    if output_file:
        write_gmt(pathways, output_file, verbose=verbose)
    return pathways
