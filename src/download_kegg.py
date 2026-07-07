"""KEGG REST API downloads and parsers."""

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
