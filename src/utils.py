"""Shared helpers: logging, caching, downloads, and output writers."""

import csv
import os
import sys
from datetime import datetime


def log_msg(*args):
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{ts}] {''.join(str(a) for a in args)}", file=sys.stderr)


def ensure_dir(path):
    os.makedirs(path, exist_ok=True)
    return path


def cached_download(url, dest, cache_dir, refresh=False):
    """Download *url* to <cache_dir>/<dest> unless it already exists."""
    ensure_dir(cache_dir)
    local = os.path.join(cache_dir, dest)
    if os.path.exists(local) and not refresh:
        log_msg("cache hit: ", dest)
        return local
    log_msg("downloading: ", url)
    import requests

    resp = requests.get(url, timeout=300)
    resp.raise_for_status()
    with open(local, "wb") as f:
        f.write(resp.content)
    if os.path.getsize(local) == 0:
        raise RuntimeError(f"Downloaded file is empty: {url}")
    return local


def write_tab(rows, columns, path):
    """Write rows (list of dicts or list of tuples) to a tab-delimited file."""
    ensure_dir(os.path.dirname(path))
    with open(path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t", quoting=csv.QUOTE_NONE,
                            escapechar="\\", quotechar=None)
        writer.writerow(columns)
        for row in rows:
            if isinstance(row, dict):
                writer.writerow([row[c] for c in columns])
            else:
                writer.writerow(row)
    log_msg("wrote: ", path, f"  ({len(rows)} rows)")


def write_pathway2gene(path2gene_pairs, out_file):
    """Write KEGG_pathway2gene.tab  (columns: v1, index)."""
    pairs = sorted(set(path2gene_pairs))
    write_tab(pairs, ["v1", "index"], out_file)


def write_pathway2name(path2name_pairs, out_file):
    """Write KEGG_pathway2name.tab  (columns: pathway, info)."""
    write_tab(path2name_pairs, ["pathway", "info"], out_file)


def write_go2gene(go2gene_pairs, out_file):
    """Write GO2gene_<NS>.tab  (columns: GO, Gene)."""
    pairs = sorted(set(go2gene_pairs))
    write_tab(pairs, ["GO", "Gene"], out_file)


def write_go2name(go2name_pairs, out_file):
    """Write GO2name_<NS>.tab  (columns: GO, Term)."""
    pairs = list(dict.fromkeys(go2name_pairs))  # deduplicate, preserve order
    write_tab(pairs, ["GO", "Term"], out_file)


def write_gmt(term2genes, term2name, out_file):
    """Write a GMT file (one line per term: TERM_ID<tab>description<tab>gene1<tab>gene2...).

    term2genes: dict mapping term_id -> list of gene_ids
    term2name:  dict mapping term_id -> term description
    """
    ensure_dir(os.path.dirname(out_file))
    count = 0
    with open(out_file, "w") as f:
        for term_id in sorted(term2genes):
            genes = sorted(set(term2genes[term_id]))
            if not genes:
                continue
            desc = term2name.get(term_id, "")
            f.write(f"{term_id}\t{desc}\t" + "\t".join(genes) + "\n")
            count += 1
    log_msg("wrote: ", out_file, f"  ({count} gene sets)")
