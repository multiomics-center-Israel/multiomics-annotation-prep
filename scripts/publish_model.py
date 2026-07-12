#!/usr/bin/env python3
"""Publish a built mummichog model as an immutable GitHub Release and record it
in the model registry (MODELS.md).

On a successful release it AUTO-APPENDS a row to MODELS.md (target/model
organism, surrogate flag, tag, asset URL, sha256, KEGG snapshot, build date) --
pulled entirely from the sidecar manifest -- and commits MODELS.md. The append
is idempotent: re-running for an already-listed tag changes nothing.

    python scripts/publish_model.py --manifest results/cre_kegg_20260711.manifest.json

The release step uses the `gh` CLI (skip it with --no-release, e.g. to only
refresh the registry). MODELS.md is tracked source, not a build artifact.
"""

import argparse
import json
import os
import shutil
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.model_registry import (  # noqa: E402
    DEFAULT_REPO,
    append_model_row,
    build_row,
    format_row,
    model_url,
    tag_from_manifest,
)
from src.utils import log_msg  # noqa: E402

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _release_notes(manifest):
    c = manifest.get("counts", {})
    v = manifest.get("validation", {})
    return (
        f"Metabolic model for `mummichog -n` "
        f"({manifest.get('model_organism', '')}, {manifest.get('source_version', '')}).\n\n"
        f"- compounds: {c.get('compounds')}, reactions: {c.get('reactions')}, "
        f"pathways: {c.get('pathways')}\n"
        f"- sha256: `{manifest.get('sha256', '')}`\n"
        f"- mummichog -n validated: {v.get('smoke_run_exit_0')}, "
        f"mass spot-check passed: {v.get('mass_spotcheck_passed')}\n"
    )


def create_release(tag, model_path, manifest_path, repo, notes):
    """Create the GitHub release + upload assets via the gh CLI.

    Tolerates an already-existing release (so the registry can still be
    refreshed). Returns True on success/already-exists, False if gh is absent.
    """
    if not shutil.which("gh"):
        log_msg("gh CLI not found; skipping release creation "
                "(use --no-release to silence, or create the release manually)")
        return False
    assets = [p for p in (model_path, manifest_path) if p and os.path.exists(p)]
    cmd = (["gh", "release", "create", tag, *assets, "--repo", repo,
            "--title", tag, "--notes", notes])
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode == 0:
        log_msg("created release: ", tag)
        return True
    if "already exists" in (proc.stderr + proc.stdout):
        log_msg("release ", tag, " already exists; continuing to registry")
        return True
    log_msg("gh release create failed: ", proc.stderr.strip())
    raise SystemExit(f"release failed for {tag}: {proc.stderr.strip()}")


def commit_models_file(models_path, tag):
    rel = os.path.relpath(models_path, REPO_ROOT)
    subprocess.run(["git", "-C", REPO_ROOT, "add", rel], check=True)
    diff = subprocess.run(["git", "-C", REPO_ROOT, "diff", "--cached", "--quiet", rel])
    if diff.returncode == 0:
        log_msg("MODELS.md unchanged; nothing to commit")
        return
    subprocess.run(["git", "-C", REPO_ROOT, "commit", "-m",
                    f"registry: add {tag} to MODELS.md"], check=True)
    log_msg("committed MODELS.md")


def main():
    parser = argparse.ArgumentParser(
        description="Publish a model release and record it in MODELS.md")
    parser.add_argument("--manifest", required=True,
                        help="Path to <stem>.manifest.json")
    parser.add_argument("--model", default=None,
                        help="Path to the model <stem>.json asset "
                             "(default: sibling of the manifest)")
    parser.add_argument("--repo", default=DEFAULT_REPO,
                        help=f"owner/repo for the release + URL (default: {DEFAULT_REPO})")
    parser.add_argument("--models-file",
                        default=os.path.join(REPO_ROOT, "MODELS.md"),
                        help="Path to MODELS.md")
    parser.add_argument("--no-release", action="store_true",
                        help="Skip the gh release; only update the registry")
    parser.add_argument("--no-commit", action="store_true",
                        help="Do not git-commit MODELS.md")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print the row that would be added; change nothing")
    args = parser.parse_args()

    with open(args.manifest) as f:
        manifest = json.load(f)
    tag = tag_from_manifest(manifest)
    model_path = args.model or os.path.join(
        os.path.dirname(os.path.abspath(args.manifest)), manifest["model_file"])

    if args.dry_run:
        log_msg("[dry-run] tag: ", tag)
        log_msg("[dry-run] url: ", model_url(manifest, args.repo))
        log_msg("[dry-run] row: ", format_row(build_row(manifest, args.repo)))
        return

    if not args.no_release:
        create_release(tag, model_path, args.manifest, args.repo,
                       _release_notes(manifest))

    added = append_model_row(args.models_file, manifest, args.repo)
    log_msg("registry: ", "added " if added else "already present, skipped ", tag)

    if added and not args.no_commit:
        commit_models_file(args.models_file, tag)


if __name__ == "__main__":
    main()
