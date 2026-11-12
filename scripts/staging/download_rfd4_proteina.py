#!/usr/bin/env python3
"""Download selected components of the nvidia/RFD4-Proteina dataset (Hugging Face), pinned to one commit.

Standard library only (no huggingface_hub), resumable, sha256-verified against the Hub's LFS digests.
Heavy: run it through sbatch (``download_rfd4_proteina.sbatch``), never on the login node.

  manifest   list the pinned revision, write MANIFEST.tsv + TASKS.tsv under the destination
  fetch      download one task (a slice of one component) from TASKS.tsv
  verify     re-check size (and, with --sha, sha256) of everything in MANIFEST.tsv

Layout under the destination mirrors the Hub repo (train_datasets/streaming-latest/...), as the
RFD4 loader expects, so the destination can be used as DATASET_STORE_ROOT.
"""

# ruff: noqa: S310  (every URL is a fixed https Hub URL)
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

REPO = "nvidia/RFD4-Proteina"
HUB = "https://huggingface.co"
PREFIX = "train_datasets/streaming-latest/"
# config name -> component directory (relative to PREFIX)
COMPONENTS = {
    "omol": "synthetic/omol_annotated",
    "rna_monomer": "synthetic/2026_09_21_rna_monomers",
    "na_complex": "synthetic/2026_09_22_protein_dna_esmfold2",
    "antibody_interfaces_distill": "synthetic/2026_09_18_antibody_antigen_esmfold2",
    "tcr_interfaces_distill": "synthetic/2026_09_09_tcr_af2",
    "synthetic_ncaa": "synthetic/ncaa/2026_09_21_ptm_esmfold2",
    "ncaa_scaffolds": "synthetic/ncaa/2026_09_26_ncaa_esmfold2_scaffolds",
    "reglyco": "synthetic/glycans/2026_09_27_reglyco",
}
ROOT_FILES = ("README.md", "dataset-inventory.json", "LICENSE")
PART_BYTES = 200 * 10**9  # split a component into tasks of about this size
CHUNK = 8 << 20


def default_dest() -> Path:
    root = os.environ.get("BIOMOL_ROOT", "/data/shared/cssb_data/BioMol_clean")
    return Path(root) / "materials" / "raw" / "rfd4_proteina"


def _get_json(url: str) -> object:
    with urllib.request.urlopen(url, timeout=120) as r:
        return json.load(r)


def cmd_manifest(a: argparse.Namespace) -> None:
    dest = Path(a.dest)
    dest.mkdir(parents=True, exist_ok=True)
    rev = a.revision or _get_json(f"{HUB}/api/datasets/{REPO}")["sha"]
    info = _get_json(f"{HUB}/api/datasets/{REPO}/revision/{rev}?blobs=true")
    files = {s["rfilename"]: s for s in info["siblings"]}
    rows, tasks = [], []
    for comp, sub in COMPONENTS.items():
        pre = f"{PREFIX}{sub}/"
        sel = sorted((p, s) for p, s in files.items() if p.startswith(pre))
        if not sel:
            sys.exit(f"no files under {pre}")
        total = sum(s["size"] for _, s in sel)
        n = max(1, -(-total // PART_BYTES))
        bins = [[0, i, []] for i in range(n)]  # greedy balance by size
        for p, s in sorted(sel, key=lambda t: -t[1]["size"]):
            b = min(bins, key=lambda x: x[0])
            b[0] += s["size"]
            b[2].append(p)
        for _, i, ps in bins:
            for p in ps:
                s = files[p]
                rows.append(
                    (comp, i, p, s["size"], (s.get("lfs") or {}).get("sha256", ""))
                )
            tasks.append((comp, i, n))
    for p in ROOT_FILES:
        if p in files:
            s = files[p]
            rows.append(
                ("root", 0, p, s["size"], (s.get("lfs") or {}).get("sha256", ""))
            )
    tasks.append(("root", 0, 1))
    (dest / "MANIFEST.tsv").write_text(
        "component\tpart\tpath\tsize\tsha256\n"
        + "".join("\t".join(map(str, r)) + "\n" for r in rows)
    )
    (dest / "TASKS.tsv").write_text(
        "task\tcomponent\tpart\tnparts\n"
        + "".join(f"{i}\t{c}\t{p}\t{n}\n" for i, (c, p, n) in enumerate(tasks))
    )
    (dest / "REVISION").write_text(rev + "\n")
    tot = sum(r[3] for r in rows)
    print(
        f"revision {rev}: {len(rows)} files, {tot / 1e9:.1f} GB, {len(tasks)} tasks -> {dest}"
    )


def read_manifest(dest: Path) -> list[tuple[str, int, str, int, str]]:
    out = []
    for line in (dest / "MANIFEST.tsv").read_text().splitlines()[1:]:
        c, part, p, size, sha = line.split("\t")
        out.append((c, int(part), p, int(size), sha))
    return out


def fetch_file(
    dest: Path, rev: str, path: str, size: int, sha: str, tries: int = 8
) -> str:
    out = dest / path
    if out.exists() and out.stat().st_size == size:
        return "skip"
    out.parent.mkdir(parents=True, exist_ok=True)
    part = out.with_name(out.name + ".part")
    url = f"{HUB}/datasets/{REPO}/resolve/{rev}/{path}"
    k = 0
    last = -1
    while k < tries:
        try:
            have = part.stat().st_size if part.exists() else 0
            if (
                have > last
            ):  # the CDN often cuts a stream short; progress is not a failure
                k, last = 0, have
            if have > size:
                part.unlink()
                have = 0
            h = hashlib.sha256()
            if have and sha:
                with open(part, "rb") as f:
                    while b := f.read(CHUNK):
                        h.update(b)
            if have < size:
                req = urllib.request.Request(
                    url, headers={"Range": f"bytes={have}-"} if have else {}
                )
                with (
                    urllib.request.urlopen(req, timeout=120) as r,
                    open(part, "ab") as f,
                ):
                    if have and r.status != 206:  # server ignored Range: start over
                        f.truncate(0)
                        have, h = 0, hashlib.sha256()
                    while b := r.read(CHUNK):
                        f.write(b)
                        h.update(b)
            if part.stat().st_size != size:
                raise OSError(f"size {part.stat().st_size} != {size}")  # noqa: TRY301  retried below
            if sha and h.hexdigest() != sha:
                part.unlink()
                raise OSError("sha256 mismatch")  # noqa: TRY301  retried below
            part.rename(out)
            return "ok"  # noqa: TRY300
        except (urllib.error.URLError, OSError, TimeoutError) as e:
            wait = min(300, 5 * 2**k)
            print(f"retry {k + 1}/{tries} {path}: {e} (sleep {wait}s)", flush=True)
            time.sleep(wait)
            k += 1
    raise RuntimeError(f"gave up on {path}")


def cmd_fetch(a: argparse.Namespace) -> None:
    dest = Path(a.dest)
    rev = (dest / "REVISION").read_text().strip()
    tasks = [ln.split("\t") for ln in (dest / "TASKS.tsv").read_text().splitlines()[1:]]
    _, comp, part, _ = tasks[a.task]
    todo = [r for r in read_manifest(dest) if r[0] == comp and r[1] == int(part)]
    total = sum(r[3] for r in todo)
    print(
        f"task {a.task}: {comp} part {part}: {len(todo)} files, {total / 1e9:.1f} GB",
        flush=True,
    )
    done = nbytes = 0
    t0 = time.time()
    with ThreadPoolExecutor(a.threads) as ex:
        futs = {
            ex.submit(fetch_file, dest, rev, p, s, sh): s for _, _, p, s, sh in todo
        }
        for f in as_completed(futs):
            f.result()
            done += 1
            nbytes += futs[f]
            if done % 100 == 0 or done == len(todo):
                print(
                    f"  {done}/{len(todo)} files, {nbytes / 1e9:.1f}/{total / 1e9:.1f} GB, "
                    f"{nbytes / 1e6 / max(1, time.time() - t0):.0f} MB/s",
                    flush=True,
                )
    print(f"task {a.task} done", flush=True)


def cmd_verify(a: argparse.Namespace) -> None:
    dest = Path(a.dest)
    rows = read_manifest(dest)
    if a.component:
        rows = [
            r
            for r in rows
            if r[0] == a.component and (a.part is None or r[1] == a.part)
        ]
    bad = 0
    for _, _, p, size, sha in rows:
        f = dest / p
        if not f.exists() or f.stat().st_size != size:
            print("MISSING/SIZE", p)
            bad += 1
        elif a.sha and sha:
            h = hashlib.sha256()
            with open(f, "rb") as fh:
                while b := fh.read(CHUNK):
                    h.update(b)
            if h.hexdigest() != sha:
                print("SHA", p)
                bad += 1
    print(
        f"verify {a.component or 'all'}{'' if a.part is None else f' part {a.part}'}: {len(rows)} files,",
        "OK" if not bad else f"{bad} bad",
    )
    sys.exit(1 if bad else 0)


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--dest", default=str(default_dest()))
    sp = ap.add_subparsers(required=True)
    m = sp.add_parser("manifest")
    m.add_argument("--revision", help="pin this commit (default: current main)")
    m.set_defaults(fn=cmd_manifest)
    f = sp.add_parser("fetch")
    f.add_argument("--task", type=int, required=True)
    f.add_argument("--threads", type=int, default=6)
    f.set_defaults(fn=cmd_fetch)
    v = sp.add_parser("verify")
    v.add_argument("--sha", action="store_true")
    v.add_argument("--component", help="only this component (default: all)")
    v.add_argument(
        "--part",
        type=int,
        help="only this part of the component (the download task's slice)",
    )
    v.set_defaults(fn=cmd_verify)
    a = ap.parse_args()
    a.fn(a)


if __name__ == "__main__":
    main()
