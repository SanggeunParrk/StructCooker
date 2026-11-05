"""Small fasta/TSV -> list projections feeding the template pipeline.

Tiny side-effect writers (like the precompute steps): each reads a fasta/TSV and
writes a derived list/fasta, returning a status string. Driven by materialize
configs whose ``output_data_path`` is the build-done marker.
"""

from __future__ import annotations

from pathlib import Path

_L_TYPE = "polypeptide(L)"


def _chain_id(header: str) -> str:
    """``100D_A_. | polypeptide(L) | ...`` -> ``100D_A`` (pdbid_chain)."""
    tok = header.split("|", 1)[0].strip()
    parts = tok.split("_")
    return "_".join(parts[:2]) if len(parts) >= 2 else tok


def filter_fasta_polypeptide_l(fasta_path: str | Path, out_path: str | Path) -> str:
    """Write only the polypeptide(L) records of ``fasta_path`` (the hmmsearch template DB)."""
    n = 0
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    keep = False
    with (
        Path(fasta_path).open(encoding="utf-8") as f,
        Path(out_path).open("w", encoding="utf-8") as out,
    ):
        for line in f:
            if line.startswith(">"):
                keep = _L_TYPE in line
                if keep:
                    out.write(line)
                    n += 1
            elif keep:
                out.write(line)
    return f"polypeptide(L) fasta: {n} records -> {out_path}"


def fasta_chain_list(fasta_path: str | Path, out_path: str | Path) -> str:
    """Write the pdbid_chain id of every polypeptide(L) record, one per line.

    This is the protein-chain work list (``template_chain_filelist``) that the
    Phase 1/2 candidate precompute and the Phase 4 keyed build iterate over. A chain
    appears once per altloc in the fasta (``1J4R_A_A``, ``1J4R_A_B``, ``1J4R_A_.``) but
    is one work item, so each id is written once, in first-seen order.
    """
    n = 0
    seen: set[str] = set()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    with (
        Path(fasta_path).open(encoding="utf-8") as f,
        Path(out_path).open("w", encoding="utf-8") as out,
    ):
        for line in f:
            if line.startswith(">") and _L_TYPE in line:
                chain = _chain_id(line[1:])
                if chain in seen:
                    continue
                seen.add(chain)
                out.write(chain + "\n")
                n += 1
    return f"chain list: {n} polypeptide(L) chains -> {out_path}"


def tsv_key_list(tsv_path: str | Path, out_path: str | Path) -> str:
    """Write column-0 keys of a TSV, one per line (e.g. seqid_to_chains -> seq_ids)."""
    n = 0
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    with (
        Path(tsv_path).open(encoding="utf-8") as f,
        Path(out_path).open("w", encoding="utf-8") as out,
    ):
        for line in f:
            key = line.split("\t", 1)[0].strip()
            if key:
                out.write(key + "\n")
                n += 1
    return f"key list: {n} keys -> {out_path}"


def glob_file_list(
    sources: list[dict], out_path: str | Path, relative_to: str | Path | None = None
) -> str:
    """Write every file matching ``{dir, glob[, recursive]}`` across ``sources``, sorted, one per line.

    For a DB whose inputs live in more than one tree -- a release's files under
    ``materials/raw`` plus files we generated under ``materials/intermediate`` -- so the
    raw tree never has to hold anything that was not downloaded. With ``relative_to`` the
    lines are relative to it (the reading config's ``data_dir``), so the list survives the
    tree being moved or renamed.
    """
    paths: list[str] = []
    for source in sources:
        root, pattern = Path(source["dir"]), str(source.get("glob", "*"))
        found = root.rglob(pattern) if source.get("recursive") else root.glob(pattern)
        paths.extend(
            str(p.relative_to(relative_to) if relative_to else p)
            for p in found
            if p.is_file()
        )
    paths.sort()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    Path(out_path).write_text("".join(f"{p}\n" for p in paths), encoding="utf-8")
    return f"file list: {len(paths)} files from {len(sources)} sources -> {out_path}"


def disordered_template_links(
    structure_dir: str | Path,
    template_dir: str | Path,
    out_dir: str | Path,
    out_list: str | Path,
    workers: int = 32,
) -> str:
    """One link per chain of the structure that has a template, and the list of the links.

    The disordered set names a template file by label chain (``templates/1a0e/1a0e_A.npz``) while its structure, and so
    ``cif_disordered``, numbers the chains (``chain_id`` 1, 2, ...). The structure also keeps one copy of each entity (1a0o:
    chains A and B), where the templates are the crystal's whole asymmetric unit (1a0o_A ... 1a0o_H, four copies of each
    entity). A consumer looks a chain up as ``<record>_<chain id>``, so the template DB must hold one key per chain id.
    For every chain of the structure the template is the file of the chain's own label if there is one, else the first
    file (label order) of the same entity; a chain with neither has no template and no key.
    ``<out_dir>/<id>/<id>_<chain id>.npz`` links to that file and the template build reads the links, keyed by their stem.
    The list holds ``<id>/<id>_<chain id>.npz`` per link. A structure that cannot be read is an error.
    """
    import os
    import shutil
    from concurrent.futures import ThreadPoolExecutor

    import numpy as np

    structure_dir, template_dir, out_dir = (
        Path(structure_dir),
        Path(template_dir),
        Path(out_dir),
    )

    def template_entity(path: Path) -> str:
        with np.load(path, allow_pickle=True) as handle:
            return str(np.unique(handle["entity_id"].astype(str))[0])

    def one(entry: str) -> tuple[list[str], dict[str, int]]:
        with np.load(
            structure_dir / entry / f"{entry}.npz", allow_pickle=True
        ) as handle:
            chains = sorted(
                set(
                    zip(
                        handle["chain_id"].astype(str).tolist(),
                        handle["label_asym_id"].astype(str).tolist(),
                        handle["entity_id"].astype(str).tolist(),
                        strict=True,
                    )
                )
            )
        files = {
            f.stem[len(entry) + 1 :]: f
            for f in sorted((template_dir / entry).glob(f"{entry}_*.npz"))
        }
        by_entity: dict[str, Path] | None = None
        links: list[str] = []
        stat = {"by_label": 0, "by_entity": 0, "unused_files": 0}
        used: set[Path] = set()
        for cid, label, entity in chains:
            src = files.get(label)
            if src is not None:
                stat["by_label"] += 1
            elif files:
                if by_entity is None:
                    by_entity = {}
                    for f in files.values():
                        by_entity.setdefault(template_entity(f), f)
                src = by_entity.get(entity)
                stat["by_entity"] += src is not None
            if src is None:
                continue
            used.add(src)
            link = out_dir / entry / f"{entry}_{cid}.npz"
            link.parent.mkdir(parents=True, exist_ok=True)
            # relative, so that renaming the release root does not break the links
            link.symlink_to(
                os.path.relpath(src, link.parent)  # relpath makes both absolute
            )
            links.append(f"{entry}/{link.name}")
        stat["unused_files"] = len(files) - len(used)
        return links, stat

    if out_dir.exists():
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True)
    entries = sorted(p.name for p in template_dir.iterdir() if p.is_dir())
    with ThreadPoolExecutor(max_workers=workers) as pool:
        results = list(pool.map(one, entries))
    links = [link for ls, _ in results for link in ls]
    keys = [Path(link).stem for link in links]
    if len(set(keys)) != len(keys):
        msg = "two template files would write the same chain key"
        raise ValueError(msg)
    totals = {
        k: sum(st[k] for _, st in results)
        for k in ("by_label", "by_entity", "unused_files")
    }
    Path(out_list).parent.mkdir(parents=True, exist_ok=True)
    Path(out_list).write_text("\n".join(links) + "\n", encoding="utf-8")
    return (
        f"disordered template links: {len(entries)} entries, {len(links)} links ({totals['by_label']} by label, "
        f"{totals['by_entity']} by entity), {totals['unused_files']} template files unused -> {out_dir}"
    )
