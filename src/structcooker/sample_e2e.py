"""Sampled end-to-end check: the whole release DAG on a few hundred inputs.

A full ``build-all`` takes days and terabytes, too much to run for every change. This
builds a small DATA_ROOT and OUTPUT_ROOT under one scratch root and runs the same 94
stages on them (``db/MANIFEST_all.yaml``, the same configs and the same code path as a
release), then compares what came out with the published DBs.

* ``prepare`` picks a deterministic sample per source -- hash-ordered ids, plus edge
  cases chosen on purpose (accepted PDB failures, manually fixed mmCIFs, repaired AFM
  heterodimers, the AFM homodimer whose MSA we generated, recovered OFD long entries, and
  PDB entries that reach train/valid/mpnn/chain) -- and mirrors their files into
  ``<root>/data`` as per-file symlinks in the production layout. Reference fastas and
  work lists are cut down to the sample; reference DBs and search databases are linked
  whole. Nothing under production is written.
* ``run`` is ``build-all`` with ``DATA_ROOT``/``OUTPUT_ROOT`` pointed at the sample (and
  every other root override cleared, so no path can resolve to production) and parallel
  searches capped at one node.
* ``verify`` reads every output and compares it with the release, key by key. Records
  that depend only on their own inputs (structures, MSAs, CCD: schemas A/E/G) must be
  byte- or value-identical; the rest (attached metadata, clusters, templates, train and
  valid sets) depend on what else is in the build, so they are reported, not required
  to match.

``prepare`` and ``verify`` scan production files and run as SLURM jobs
(``structcooker e2e-sample``); see docs/e2e-build.md.
"""
from __future__ import annotations

import hashlib
import heapq
import json
import os
import re
import sys
from collections import Counter
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator

PROD_DATA_ROOT = Path("/data/shared/cssb_data")
PROD_OUTPUT_ROOT = PROD_DATA_ROOT / "BioMol_clean"
MATERIALS = Path("BioMol/materials")
SEED = "structcooker-e2e-sample-v1"

# Env vars the db configs consult before DATA_ROOT/OUTPUT_ROOT. A sampled run clears them
# all so every path resolves under the sample root.
ROOT_OVERRIDES = (
    "MMCIF_ROOT", "DISTILLATION_ROOT", "OFD_INTERMEDIATE_ROOT", "TEDDYMER_ROOT", "AFDB_MULTIMER_ROOT",
    "MSA_ROOT", "REFERENCE_ROOT", "SEQ_CLUSTER30_PATH", "MONOMER_SEQID_MAP", "DISORDERED_SEQID_MAP",
    "OPENFOLD_STRUCTURE_OVERRIDES",
)
DEFAULT_COUNTS = {"pdb": 30, "ofd_short": 20, "ofd_long": 20, "ofd_rna": 20, "ofd_disordered": 20,
                  "ofd_disordered_templates": 20, "teddymer": 20, "afm_homodimer": 20, "afm_heterodimer": 20}
# Random PDB picks above this size are skipped (a capsid's attached record is tens of GB).
PDB_MAX_BYTES = 2_000_000
# Published PDB-derived DBs: a couple of their entries are put in the sample so every
# downstream filter (train cutoffs, valid splits, mpnn, chain) has something to keep.
PDB_SEED_DBS = (
    "pdb/cif_attached/cif_pdb_attached_train_20210930_res9_chain300.lmdb",
    "pdb/cif_attached/cif_pdb_attached_train_20260301_res9_chain300.lmdb",
    "pdb/cif/cif_mpnn_train_20210930_res9_chain300.lmdb",
    "pdb/valid/cif_valid_1.lmdb",
    "pdb/valid/cif_attached_valid_2.lmdb",
    "chain/cif/cif_chain.lmdb",
)
PDB_SEEDS_PER_DB = 2
MANUAL_FIXES_IN_SAMPLE = 2
# Two of the 14 OFD long entries recovered on 2026-09-14 (docs/long-cif-recovery-2026-09-14.md):
# absent from the reference long.fasta, so they exercise ofd_long_recovered_fasta.
OFD_LONG_RECOVERED = ("MGYP000485802788", "MGYP006508416081")
# The AFM homodimer whose MSA the release lacks; we generate it (docs/e2e-build.md §2).
AFM_GENERATED_MSA_ENTITY = "AF-0000000208803231"
# Schemas whose record is a function of that record's own inputs: must equal the release.
EXACT_SCHEMAS = frozenset({"A", "E", "G"})


# --------------------------------------------------------------------------- helpers
def pick(items: Iterable[str], n: int, salt: str = "") -> list[str]:
    """The ``n`` items first in seeded-hash order: deterministic, spread over the set.

    Picks from overlapping sets need different salts: unsalted, the first train entries
    are the first PDB entries too, and a seed pick would add nothing to the random one.
    """
    def rank(item: str) -> bytes:
        return hashlib.blake2b(f"{SEED}:{salt}:{item}".encode(), digest_size=8).digest()
    return sorted(heapq.nsmallest(n, set(items), key=rank))


def link(src: Path, dst: Path) -> bool:
    """Symlink one file into the sample tree; False when the source does not exist."""
    if not src.is_file():
        return False
    dst.parent.mkdir(parents=True, exist_ok=True)
    if not dst.is_symlink():
        dst.symlink_to(src)
    return True


def link_tree(src: Path, dst: Path, pattern: str = "*") -> int:
    """Mirror a directory as real directories holding file symlinks.

    Input fingerprints refuse directory symlinks below an input root, so only files are
    linked.
    """
    n = 0
    for path in sorted(src.rglob(pattern)) if src.is_dir() else []:
        if path.is_file():
            n += link(path, dst / path.relative_to(src))
    return n


def read_fasta(path: Path) -> Iterator[tuple[str, str]]:
    """``(header without '>', sequence)`` pairs."""
    header, chunks = None, []
    with path.open() as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header, chunks = line[1:], []
            elif line:
                chunks.append(line)
    if header is not None:
        yield header, "".join(chunks)


def fasta_entry_id(header: str) -> str:
    """``MGYP..._1_. | ...`` / ``AAAA02045115.1_275-1799_1_. | ...`` -> the entry id."""
    return header.split()[0].rsplit("_", 2)[0]


def write_fasta_subset(src: Path, dst: Path, ids: set[str]) -> dict[str, str]:
    """Copy the records of ``ids`` from ``src``; returns header -> sequence written."""
    kept = {h: s for h, s in read_fasta(src) if fasta_entry_id(h) in ids}
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text("".join(f">{h}\n{s}\n" for h, s in kept.items()))
    return kept


def index_keys(index: Path) -> Iterator[str]:
    """Keys of a published ``<db>.index.tsv``."""
    with index.open() as handle:
        for line in handle:
            if line.startswith("#") or line.startswith("key\t"):
                continue
            yield line.split("\t", 1)[0]


def tsv_rows(path: Path) -> Iterator[list[str]]:
    with path.open() as handle:
        for line in handle:
            yield line.rstrip("\n").split("\t")


# --------------------------------------------------------------------------- prepare
class Sampler:
    """Builds ``<root>/data`` (DATA_ROOT) and ``<root>/out`` (OUTPUT_ROOT) for a sample."""

    def __init__(self, root: Path, counts: dict[str, int],
                 prod_data: Path = PROD_DATA_ROOT, prod_out: Path = PROD_OUTPUT_ROOT) -> None:
        self.root, self.counts = root, counts
        self.data, self.out = root / "data", root / "out"
        self.prod_data, self.prod_out = prod_data, prod_out
        self.raw, self.inter = prod_data / MATERIALS / "raw", prod_data / MATERIALS / "intermediate"
        self.s_raw, self.s_inter = self.data / MATERIALS / "raw", self.data / MATERIALS / "intermediate"
        self.report: dict[str, Any] = {}

    def run(self) -> dict[str, Any]:
        if (self.root / "sample.json").exists():
            msg = f"{self.root} already holds a sample; use a new root"
            raise FileExistsError(msg)
        self.shared()
        self.pdb()
        self.ofd()
        self.teddymer()
        self.afm()
        self.report["counts"] = self.counts
        (self.root / "sample.json").write_text(json.dumps(self.report, indent=1, sort_keys=True))
        return self.report

    def shared(self) -> None:
        """Inputs linked whole: search DBs, CCD, the seq_id seed and reference DBs."""
        (self.data / "seq_db").parent.mkdir(parents=True, exist_ok=True)
        (self.data / "seq_db").symlink_to(self.prod_data / "seq_db")
        link(self.prod_data / "external/SabDab/sabdab_summary_all.tsv",
             self.data / "external/SabDab/sabdab_summary_all.tsv")
        link((self.prod_data / "reference/seq_id_map.tsv").resolve(), self.data / "reference/seq_id_map.tsv")
        link_tree(self.inter / "ccd/components", self.s_inter / "ccd/components")
        link_tree(self.inter / "openfold_distillation/long_recovered",
                  self.s_inter / "openfold_distillation/long_recovered")
        (self.out / "metadata").mkdir(parents=True, exist_ok=True)
        for name in ("reference_ccd.lmdb", "reference_chain_selection.lmdb"):
            link_tree(self.prod_out / "metadata" / name, self.out / "metadata" / name)
        (self.out / "materials/intermediate/mmseqs_tmp").mkdir(parents=True, exist_ok=True)

    # -- PDB --------------------------------------------------------------------------
    def pdb(self) -> None:
        cif_root = self.raw / "cif"
        sizes = {p.name.split(".")[0]: p.stat().st_size
                 for d in cif_root.iterdir() if d.is_dir() for p in d.iterdir() if p.name.endswith(".cif.gz")}
        small = [i for i, s in sizes.items() if s <= PDB_MAX_BYTES]
        chosen = {i: "random" for i in pick(small, self.counts["pdb"])}
        for db in PDB_SEED_DBS:
            ids = {k[:4] for k in index_keys(self.prod_out / "lmdb" / f"{db}.index.tsv")}
            for i in pick((i for i in ids if sizes.get(i, PDB_MAX_BYTES + 1) <= PDB_MAX_BYTES), PDB_SEEDS_PER_DB, salt=db):
                chosen.setdefault(i, f"seed:{db}")
        policy = json.loads((Path(__file__).resolve().parents[2] / "db/pdb-exceptions.json").read_text())
        for i in policy.get("pdb/cif", {}):
            chosen.setdefault(i, "accepted failure")
        fixes = Path(__file__).resolve().parents[2] / "db/pdb/manual_cif_fixes.txt"
        fixed = [s.strip() for s in fixes.read_text().splitlines() if s.strip() and not s.startswith("#")]
        for i in fixed[:MANUAL_FIXES_IN_SAMPLE]:
            chosen.setdefault(i, "manual cif fix")
        for i in chosen:
            link(cif_root / i[1:3] / f"{i}.cif.gz", self.s_raw / "cif" / i[1:3] / f"{i}.cif.gz")

        # Their chains' seq_ids (published fasta + map) -> the provided MSAs and SignalP files.
        upper = {i.upper() for i in chosen}
        seqs = {s for h, s in read_fasta(self.prod_out / "materials/intermediate/fasta/cif_pdb.fasta")
                if h.split("_", 1)[0] in upper}
        seq_ids = [sid for sid, seq, *_ in tsv_rows(self.prod_out / "metadata/seq_id_map.tsv") if seq in seqs]
        msas = 0
        for sid in seq_ids:
            if sid[0] in "PR":
                msas += link_tree(self.inter / "msa" / sid[0] / sid[-3:] / sid,
                                  self.s_inter / "msa" / sid[0] / sid[-3:] / sid) > 0
                link(self.prod_out / f"materials/intermediate/signalp/{sid}.gff3",
                     self.out / f"materials/intermediate/signalp/{sid}.gff3")
        (self.out / "materials/intermediate/signalp").mkdir(parents=True, exist_ok=True)
        self.report["pdb"] = {"entries": dict(sorted(chosen.items())), "seq_ids": len(seq_ids), "msas": msas}

    # -- OpenFold distillation --------------------------------------------------------
    def ofd(self) -> None:
        ofd_raw = self.raw / "openfold_distillation"
        s_ofd_raw = self.s_raw / "openfold_distillation"
        fasta, s_fasta = self.inter / "openfold_distillation/fasta", self.s_inter / "openfold_distillation/fasta"
        for subset in ("short", "long"):
            ids = set(pick((fasta_entry_id(h) for h, _ in read_fasta(fasta / f"{subset}.fasta")),
                           self.counts[f"ofd_{subset}"]))
            write_fasta_subset(fasta / f"{subset}.fasta", s_fasta / f"{subset}.fasta", ids)
            extra = set(OFD_LONG_RECOVERED) if subset == "long" else set()
            base = Path("monomer_distillation_sets_v2") / f"{subset}_monomers"
            models = 0
            for i in sorted(ids | extra):
                for name in ("best_structure_relaxed.pdb", "best_structure_relaxed.pdb.zst"):
                    models += link(ofd_raw / base / "raw" / i / name, s_ofd_raw / base / "raw" / i / name)
                link_tree(ofd_raw / base / "preprocessed" / i, s_ofd_raw / base / "preprocessed" / i)
            self.report[f"ofd_{subset}"] = {"entries": sorted(ids), "recovered": sorted(extra), "models": models}

        rna = set(pick((fasta_entry_id(h) for h, _ in read_fasta(fasta / "rna.fasta")), self.counts["ofd_rna"]))
        write_fasta_subset(fasta / "rna.fasta", s_fasta / "rna.fasta", rna)
        cache = Path("rna_distillation_set/rna_monomer_preprocessed_cache")
        for i in rna:
            link_tree(ofd_raw / cache / i, s_ofd_raw / cache / i)
        self.report["ofd_rna"] = {"entries": sorted(rna)}

        dis = set(pick((fasta_entry_id(h) for h, _ in read_fasta(fasta / "disordered.fasta")),
                       self.counts["ofd_disordered"]))
        write_fasta_subset(fasta / "disordered.fasta", s_fasta / "disordered.fasta", dis)
        for i in dis:
            link_tree(ofd_raw / "disordered_set/structure_files" / i, s_ofd_raw / "disordered_set/structure_files" / i)
        work, s_work = self.inter / "openfold_distillation/disordered_set", self.s_inter / "openfold_distillation/disordered_set"
        s_work.mkdir(parents=True, exist_ok=True)

        def ours(stem: str) -> bool:
            return stem.rsplit("_", 1)[0] in dis
        msa_lines = [s for s in (work / "msa_filelist.txt").read_text().splitlines() if s and ours(Path(s).stem)]
        (s_work / "msa_filelist.txt").write_text("".join(f"{s}\n" for s in msa_lines))
        (s_work / "msa_stem_seqid.tsv").write_text(
            "".join("\t".join(r) + "\n" for r in tsv_rows(work / "msa_stem_seqid.tsv") if ours(r[0])))
        templates = pick((s for s in (work / "protein_tmpl_filelist.txt").read_text().splitlines() if s),
                         self.counts["ofd_disordered_templates"])
        (s_work / "protein_tmpl_filelist.txt").write_text("".join(f"{s}\n" for s in templates))
        self.report["ofd_disordered"] = {"entries": sorted(dis), "msas": len(msa_lines), "templates": len(templates)}

    # -- teddymer ---------------------------------------------------------------------
    def teddymer(self) -> None:
        root = self.raw / "teddymer"
        names = {r[0]: Path(r[1]).name for r in tsv_rows(root / "SOURCE.tsv") if r[0] != "dimer_index"}
        chosen = pick(names, self.counts["teddymer"])
        for i in chosen:
            link(root / "pdb" / i[-3:] / names[i], self.s_raw / "teddymer/pdb" / i[-3:] / names[i])
        self.report["teddymer"] = {"entries": [names[i] for i in chosen]}

    # -- AFDB multimer ----------------------------------------------------------------
    def afm(self) -> None:
        root, s_root = self.raw / "afdb_multimer", self.s_raw / "afdb_multimer"
        afm_inter, s_afm_inter = self.inter / "afdb_multimer", self.s_inter / "afdb_multimer"

        def entities(source: Path) -> list[str]:
            return [r[0] for r in tsv_rows(source) if r[0] != "entity_id"]

        def link_msa(entity: str) -> bool:
            name = f"{entity}-msa_v1.a3m.zst"
            return link(root / "msa" / entity[-3:] / name, s_root / "msa" / entity[-3:] / name)

        homo = sorted({*pick(entities(root / "homodimer/SOURCE.tsv"), self.counts["afm_homodimer"]),
                       AFM_GENERATED_MSA_ENTITY})
        for e in homo:
            link(root / "homodimer/pdb" / e[-3:] / f"{e}-model_v1.pdb", s_root / "homodimer/pdb" / e[-3:] / f"{e}-model_v1.pdb")
            link_msa(e)

        overrides = json.loads((afm_inter / "heterodimer/repaired/overrides.json").read_text())
        repaired = {Path(p).name.removesuffix("-model_v1.cif.gz"): p for p in overrides}
        chain_msa = list(tsv_rows(afm_inter / "heterodimer/chain_msa.tsv"))
        header, rows = chain_msa[0], chain_msa[1:]
        hetero = sorted(set(pick(entities(root / "heterodimer/SOURCE.tsv"), self.counts["afm_heterodimer"]))
                        | set(repaired))
        sample_overrides = {}
        for e in hetero:
            src = root / "heterodimer/cif" / e[-3:] / f"{e}-model_v1.cif.gz"
            dst = s_root / "heterodimer/cif" / e[-3:] / f"{e}-model_v1.cif.gz"
            link(src, dst)
            if e in repaired:
                sample_overrides[str(dst)] = overrides[repaired[e]]
        kept = [r for r in rows if r[0] in set(hetero)]
        (s_afm_inter / "heterodimer/repaired").mkdir(parents=True, exist_ok=True)
        (s_afm_inter / "heterodimer/chain_msa.tsv").write_text("".join("\t".join(r) + "\n" for r in [header, *kept]))
        (s_afm_inter / "heterodimer/repaired/overrides.json").write_text(json.dumps(sample_overrides, indent=1))
        linked = sum(link_msa(r[3]) for r in kept)
        link_tree(afm_inter / "msa_generated", s_afm_inter / "msa_generated", "AF-*")
        self.report["afm"] = {"homodimer": homo, "heterodimer": hetero, "repaired": sorted(set(repaired) & set(hetero)),
                              "heterodimer_msas": linked}


# --------------------------------------------------------------------------- run
def sample_environment(root: Path) -> dict[str, str]:
    """``os.environ`` for a sampled run: roots at the sample, every override cleared."""
    env = {k: v for k, v in os.environ.items() if k not in ROOT_OVERRIDES}
    env.update(DATA_ROOT=str(root / "data"), OUTPUT_ROOT=str(root / "out"), STRUCTCOOKER_MAX_NODES="1")
    # Nodes whose /dev/shm cannot hold joblib's semaphores (docs/e2e-build.md §4).
    env.setdefault("DATACOOKER_EXCLUDE_NODES", "node01,node04,node07")
    return env


def escaping_paths(configs: dict[str, dict], root: Path) -> list[str]:
    """Inputs or outputs of a sampled run that resolve outside the sample root."""
    from structcooker import preflight

    bad = []
    for name, cfg in configs.items():
        for path in [*preflight.input_paths(cfg), *preflight.output_paths(cfg)]:
            if not Path(path).is_relative_to(root):
                bad.append(f"{name}: {path}")
    return bad


# --------------------------------------------------------------------------- verify
def diff(a: Any, b: Any, path: str = "") -> str | None:
    """First difference between two decoded records, as a path; None when equal."""
    import numpy as np

    if isinstance(a, dict) and isinstance(b, dict):
        if a.keys() != b.keys():
            return f"{path}{{keys}}"
        for k in a:
            found = diff(a[k], b[k], f"{path}/{k}")
            if found:
                return found
        return None
    if isinstance(a, np.ndarray) or isinstance(b, np.ndarray):
        a, b = np.asarray(a), np.asarray(b)
        if a.shape != b.shape:
            return f"{path}[shape]"
        if a.dtype.kind == "f":
            return None if np.array_equal(a, b, equal_nan=True) else path
        return None if np.array_equal(a, b) else path
    if isinstance(a, (list, tuple)) and isinstance(b, (list, tuple)):
        if len(a) != len(b):
            return f"{path}[len]"
        for i, (x, y) in enumerate(zip(a, b, strict=True)):
            found = diff(x, y, f"{path}/{i}")
            if found:
                return found
        return None
    if hasattr(a, "__dict__") and type(a) is type(b):
        return diff(vars(a), vars(b), f"{path}<{type(a).__name__}>")
    try:
        return None if a == b else path
    except ValueError:                           # array-likes without a truth value
        return None if repr(a) == repr(b) else path


def seq_id_translation(sample_map: Path, release_map: Path) -> dict[str, str]:
    """Sample seq_id -> release seq_id, by sequence, where the two differ.

    Sequences missing from the seq_id seed (teddymer, AFM, recovered OFD long) are numbered
    in build order, so a sample numbers them differently from the release.
    """
    # Keyed by type prefix too: a protein and an RNA can share a string ("ACGU...").
    by_seq = {(sid[0], seq): sid for sid, seq, *_ in tsv_rows(sample_map)}
    out = {}
    for sid, seq, *_ in tsv_rows(release_map):
        mine = by_seq.get((sid[0], seq))
        if mine is not None and mine != sid:
            out[mine] = sid
    return out


def translate(value: Any, table: dict[str, str]) -> Any:
    """``value`` with every string that is a sample seq_id replaced by the release one."""
    import numpy as np

    if not table:
        return value
    if isinstance(value, dict):
        return {translate(k, table): translate(v, table) for k, v in value.items()}
    if isinstance(value, list):
        return [translate(v, table) for v in value]
    if isinstance(value, tuple):
        return tuple(translate(v, table) for v in value)
    if isinstance(value, str):
        return table.get(value, value)
    if isinstance(value, np.ndarray) and value.dtype.kind in "UO":
        return np.vectorize(lambda v: table.get(v, v) if isinstance(v, str) else v, otypes=[value.dtype])(value) \
            if value.size else value
    return value


def compare_lmdb(sample: Path, prod: Path, schema: str, table: dict[str, str] | None = None) -> dict[str, Any]:
    """Every sample record against the release record of the same key.

    ``table`` (:func:`seq_id_translation`) maps the sample's own seq_ids to the release's,
    in keys and in record values.
    """
    table = table or {}
    from datacooker.lmdb.sharded import open_env

    from structcooker import schemas
    from structcooker.instructions.transforms.codecs import from_bytes

    raw = schemas.get(schema).codec == schemas.Codec.raw
    counts: Counter[str] = Counter()
    paths: Counter[str] = Counter()
    examples: dict[str, str] = {}
    kw = {"readonly": True, "lock": False, "readahead": False}
    with open_env(str(sample), **kw) as s_env, open_env(str(prod), **kw) as p_env, \
            s_env.begin() as s_txn, p_env.begin() as p_txn:
        for key, value in s_txn.cursor():
            counts["records"] += 1
            key = table.get(key.decode(), key.decode()).encode()
            other = p_txn.get(key)
            if other is None:
                counts["absent_in_release"] += 1
                examples.setdefault("absent_in_release", key.decode())
                continue
            if bytes(value) == bytes(other):
                counts["identical"] += 1
                continue
            found = None if raw else diff(translate(from_bytes(bytes(value)), table), from_bytes(bytes(other)))
            if found is None and not raw:
                counts["identical"] += 1          # same values, different encoding
                continue
            counts["different"] += 1
            pattern = re.sub(r"/\d+", "/#", found or "<bytes>")
            paths[pattern] += 1
            examples.setdefault(pattern, key.decode())
    return {**counts, "differing_fields": dict(paths.most_common(10)), "examples": examples}


def compare_lines(sample: Path, prod: Path, table: dict[str, str] | None = None) -> dict[str, Any]:
    """How many lines of a sample projection appear verbatim in the release's.

    Renumbered seq_ids (``table``) are translated first, field by field.
    """
    table = table or {}
    lines = {"\t".join(table.get(f, f) for f in line.split("\t")) if table else line
             for line in sample.read_text().splitlines()}
    found = 0
    with prod.open() as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line in lines:
                found += 1
                lines.discard(line)
    return {"lines": found + len(lines), "in_release": found}


def verify(root: Path, run_dir: Path, manifest: Path) -> dict[str, Any]:
    """Per node: did the stage complete, and how do its outputs compare with the release."""
    from datacooker.config import load_config

    from structcooker import preflight, workflow_ops

    repo = Path(__file__).resolve().parents[2]
    env = sample_environment(root)
    os.environ.clear()
    os.environ.update(env)
    import yaml

    graph = yaml.safe_load(manifest.read_text())
    summary = json.loads((run_dir / "summary.json").read_text()) if (run_dir / "summary.json").exists() else {}
    states = summary.get("nodes", {})
    out_root = root / "out"
    sample_map = out_root / "metadata/seq_id_map.tsv"
    table = seq_id_translation(sample_map, PROD_OUTPUT_ROOT / "metadata/seq_id_map.tsv") if sample_map.exists() else {}
    nodes: dict[str, Any] = {"_seq_id_translation": {"renumbered": len(table)}}
    failures: list[str] = []
    for name in graph:
        cfg = load_config(repo / "db" / f"{name}.yaml")
        state = json.loads((run_dir / "nodes" / name / "state.json").read_text()) \
            if (run_dir / "nodes" / name / "state.json").exists() else {}
        entry: dict[str, Any] = {"phase": states.get(name, state.get("phase")), "op": workflow_ops.op_of(cfg)}
        if state.get("error"):
            entry["error"] = state["error"]
        if entry["phase"] != "complete":
            failures.append(f"{name}: {entry['phase']}")
        target = cfg.get("new_env_path") or cfg.get("env_path")
        try:
            if target and Path(target).exists():
                schema = str(cfg["schema"])
                prod = PROD_OUTPUT_ROOT / Path(target).relative_to(out_root)
                entry["schema"] = schema
                entry["exact"] = schema in EXACT_SCHEMAS
                entry.update(compare_lmdb(Path(target), prod, schema, table) if prod.exists() else {"release": "absent"})
                if entry["exact"] and (entry.get("different") or entry.get("absent_in_release")):
                    failures.append(f"{name}: {entry.get('different', 0)} different, "
                                    f"{entry.get('absent_in_release', 0)} absent from the release")
            elif not target:
                files = {}
                for out in preflight.output_paths(cfg):
                    prod = PROD_OUTPUT_ROOT / out.relative_to(out_root)
                    if out.is_file() and prod.is_file():
                        files[str(out.relative_to(out_root))] = compare_lines(out, prod, table)
                    elif out.is_dir():
                        files[str(out.relative_to(out_root))] = {"files": sum(1 for p in out.rglob("*") if p.is_file())}
                entry["outputs"] = files
        except Exception as exc:                   # noqa: BLE001 - one node's report must not stop the rest
            entry["verify_error"] = f"{type(exc).__name__}: {exc}"
            failures.append(f"{name}: verify error {exc}")
        nodes[name] = entry
    report = {"complete": bool(summary.get("complete")), "failures": failures, "nodes": nodes}
    (root / "verify.json").write_text(json.dumps(report, indent=1, sort_keys=True))
    return report


def main() -> None:
    """``python -m structcooker.sample_e2e prepare|verify ...`` -- inside a SLURM job."""
    if not os.environ.get("SLURM_JOB_ID"):
        msg = "Sampled e2e prepare/verify scan production files; submit them with `structcooker e2e-sample`"
        raise RuntimeError(msg)
    action, root = sys.argv[1], Path(sys.argv[2])
    if action == "prepare":
        counts = json.loads(sys.argv[3]) if len(sys.argv) > 3 else DEFAULT_COUNTS
        report = Sampler(root, {**DEFAULT_COUNTS, **counts}).run()
        print(json.dumps({k: (len(v) if isinstance(v, (list, dict)) else v) for k, v in report.items()}))
    elif action == "verify":
        report = verify(root, Path(sys.argv[3]), Path(sys.argv[4]))
        print(json.dumps({"complete": report["complete"], "failures": report["failures"]}, indent=1))
        if report["failures"]:
            sys.exit(1)
    else:
        msg = f"unknown action {action}"
        raise SystemExit(msg)


if __name__ == "__main__":
    main()
