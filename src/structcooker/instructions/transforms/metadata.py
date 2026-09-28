import logging
import os
import tempfile
import time
from pathlib import Path
from typing import cast

import numpy as np
from biomol.core import NodeFeature

from structcooker.instructions.transforms.sequence import filter_water
from structcooker.mols import CIFMol, CIFMolAttached
from structcooker.utils.mapping import cluster_id, mol_type_map

_SEQMETA_SEP = "\t"  # separates seq_id and cluster in an LmdbDict value
logger = logging.getLogger(__name__)


class LmdbDict:
    """Read-only ``cif_id -> (seq_id, seq_cluster)`` map backed by an on-disk LMDB.

    The seq metadata map over ALL PDB chains is multi-GB as a Python dict. Shipping it
    inside every Ray worker's state made each of the ~112 workers on a node deserialize
    its own multi-GB copy -> hundreds of GB -> OOM (even on tiny structures, since the
    footprint was per-worker metadata, not per-item). This wraps the map as an LMDB
    instead: the shard driver writes it once, and every worker on the node mmaps the
    SAME file (shared OS page cache), so per-worker resident metadata is ~0. It pickles
    to just its path, so shipping it through Ray costs nothing. Behaves like the dict it
    replaces for the two operations attach_metadata needs: ``in`` and ``[]``.
    """

    def __init__(self, path: str | Path) -> None:
        self._path = str(path)
        self._env = None  # opened lazily per worker process, never pickled

    def __getstate__(self) -> dict:
        """Pickle only the path -- the env handle is per-process, reopened lazily."""
        return {"_path": self._path}

    def __setstate__(self, state: dict) -> None:
        """Restore from a pickle: keep the path, defer opening the env."""
        self._path = state["_path"]
        self._env = None

    def _get_env(self):  # noqa: ANN202
        if self._env is None:
            from datacooker.lmdb.sharded import open_env

            self._env = open_env(
                self._path, readonly=True, lock=False, subdir=True, max_readers=4096,
            )
        return self._env

    def __contains__(self, key: str) -> bool:
        """Return True if ``key`` (a cif_id) has metadata in the LMDB."""
        with self._get_env().begin() as txn:
            return txn.get(key.encode()) is not None

    def __getitem__(self, key: str) -> tuple[str, str]:
        """Return ``(seq_id, seq_cluster)`` for ``key``; raise KeyError if absent."""
        with self._get_env().begin() as txn:
            raw = txn.get(key.encode())
        if raw is None:
            raise KeyError(key)
        seqid, cluster = raw.decode().split(_SEQMETA_SEP)
        return seqid, cluster


def seq_metadata_map_to_lmdb(seq_metadata_map: dict[str, tuple[str, str]]) -> LmdbDict:
    """Write the seq-metadata dict to a node-local LMDB and return an ``LmdbDict``.

    Runs once on the shard driver (where the dict is built); every Ray worker on the
    node then mmaps this file instead of holding its own copy. The file lands under
    ``TMPDIR`` (node-local scratch) so each shard's workers share it via page cache.
    """
    import lmdb

    scratch = tempfile.mkdtemp(prefix="seqmeta_", dir=os.environ.get("TMPDIR", "/tmp"))  # noqa: S108
    path = os.path.join(scratch, "seq_metadata.lmdb")  # noqa: PTH118
    env = lmdb.open(path, map_size=32 * (1 << 30), subdir=True)
    with env.begin(write=True) as txn:
        for cif_id, (seqid, cluster) in seq_metadata_map.items():
            txn.put(cif_id.encode(), f"{seqid}{_SEQMETA_SEP}{cluster}".encode())
    env.close()
    return LmdbDict(path)


def load_tsv(
    tsv_file_path: Path,
    *,
    split_by_comma: bool = True,
) -> dict[str, list[str]]:
    """Load a TSV file and return its contents as a dictionary."""
    data: dict[str, list[str]] = {}
    with tsv_file_path.open("r", encoding="utf-8") as file:
        lines = file.readlines()
        for line in lines:
            key, value = line.strip().split("\t")
            if not split_by_comma:
                data[key] = [value]
            else:
                data[key] = value.split(",")
    return data


def reverse_dict(
    input_dict: dict[str, list[str]],
) -> dict[str, str]:
    """Reverse a dictionary mapping from str to list[str] into a dictionary mapping from str to str."""
    output_dict: dict[str, str] = {}
    for key, value_list in input_dict.items():
        for value in value_list:
            output_dict[value] = key
    return output_dict


def load_pairs(tsv_file_path: Path) -> dict[str, list[str]]:
    """Read a two-column relation without dropping repeated left-hand keys."""
    pairs: dict[str, list[str]] = {}
    with tsv_file_path.open(encoding="utf-8") as source:
        for line in source:
            left, right = line.rstrip("\r\n").split("\t")
            pairs.setdefault(left, []).append(right)
    return pairs


def build_seqid_map(
    seqid2seq: dict[str, list[str]],
) -> dict[str, dict[str, str]]:
    """Build a mapping from sequence and molecule type to sequence ID."""
    seqid_map: dict[str, dict[str, str]] = {}  # first key : mol type, second key : seq
    for mol_identifier in mol_type_map.values():
        seqid_map[mol_identifier] = {}
    for seqid, seq_list in seqid2seq.items():
        seq = seq_list[0]
        mol_identifier = seqid[0]
        if mol_identifier not in seqid_map:
            msg = f"Unknown molecule identifier {mol_identifier} in seqid {seqid}."
            raise ValueError(msg)
        seqid_map[mol_identifier][seq] = seqid
    return seqid_map


def build_seq_metadata_map(
    raw_fasta_dict: dict[str, str],
    seqid2seq: dict[str, list[str]],
    seqclusters2seqids: dict[str, list[str]],
    db_code: str,
) -> dict[str, tuple[str, str]]:
    """Build a metadata map from sequence cluster ID to sequence ID."""
    seq_metadata_map: dict[str, tuple[str, str]] = {}  # cif_id -> (seq id, seq cluster)
    seqid2seqcluster: dict[str, str] = reverse_dict(seqclusters2seqids)
    seqid_map = build_seqid_map(seqid2seq)

    for header, sequence in raw_fasta_dict.items():
        mol_type = header.split("|")[1].strip()
        mol_identifier = mol_type_map.get(mol_type, "X")
        seqid = seqid_map[mol_identifier].get(sequence)
        if seqid is None:
            msg = f"Sequence ID not found for molecule type {mol_identifier} and sequence {sequence}."
            raise KeyError(msg)
        # A seq_id absent from this DB's clustering is a sequence that clustering did not
        # contain -- e.g. a from-scratch id past the seeded universe -- so it forms its own
        # singleton cluster (the seq_id is its own representative), rather than erroring.
        seqcluster = seqid2seqcluster.get(seqid, cluster_id(db_code, seqid))
        cif_id = header.split("|")[0].strip()  # pdbid_chainid_altid
        seq_metadata_map[cif_id] = (seqid, seqcluster)
    return seq_metadata_map


def parse_signalp(signalp_path: Path) -> tuple[int, int] | None:
    """Parse the signalp output file and extract the sequence ids."""
    if not signalp_path.exists():
        return None
    with signalp_path.open("r", encoding="utf-8") as f:
        lines = f.readlines()

    result = lines[1].split("\t")
    return int(result[3]) - 1, int(result[4]) - 1


def load_signalp(
    signalp_dir: Path | None,
) -> dict[str, tuple[int, int]]:
    """Load SignalP data from a directory containing GFF3 files."""
    signalp_data = {}
    if signalp_dir is not None:
        if not signalp_dir.exists():
            msg = f"SignalP directory {signalp_dir} does not exist."
            raise FileNotFoundError(msg)
        for signalp_file in signalp_dir.glob("*.gff3"):
            seqid = signalp_file.stem
            signalp_data[seqid] = parse_signalp(signalp_file)
    return signalp_data


def index_signalp_sequences(
    seqid2seq_path: Path,
    signalp_dict: dict[str, tuple[int, int]],
) -> dict[tuple[str, str], tuple[int, int]]:
    """Stream the sequence table once, retaining only signal-peptide predictions."""
    indexed: dict[tuple[str, str], tuple[int, int]] = {}
    if not signalp_dict:
        return indexed
    predictions: dict[str, tuple[int, int]] = {}
    for seqid, prediction in signalp_dict.items():
        key = _sequence_counter_key(seqid)
        if key in predictions and predictions[key] != prediction:
            msg = f"Conflicting SignalP predictions for sequence counter {key}."
            raise ValueError(msg)
        predictions[key] = prediction
    with seqid2seq_path.open(encoding="utf-8") as source:
        for line in source:
            seqid, sequence = line.rstrip("\r\n").split("\t")
            prediction = predictions.get(_sequence_counter_key(seqid))
            if prediction is not None:
                indexed[seqid[0], sequence] = prediction
            else:
                # Match build_seqid_map: the final ID for a sequence wins.
                indexed.pop((seqid[0], sequence), None)
    if not indexed:
        logger.warning(
            "No SignalP prediction IDs matched %s; check prediction provenance and ID format.",
            seqid2seq_path,
        )
    return indexed


def _sequence_counter_key(seqid: str) -> str:
    """Keep molecule identity while ignoring zero padding on assigned counters."""
    if len(seqid) > 1 and seqid[1:].isdecimal():
        return seqid[0] + str(int(seqid[1:]))
    return seqid


def attach_metadata(
    cifmol: CIFMol,
    seq_metadata_map: "dict[str, tuple[str, str]] | LmdbDict",
) -> dict | None:
    """Attach metadata to a CIFMol object.

    Water is removed first so the chain set matches the fasta / seq_id_map
    (``build_fasta`` filters water before assigning sequence ids); an all-water
    entry yields ``None`` and is skipped by the rebuild runner.
    """
    cifmol_nw = filter_water(cifmol)
    if cifmol_nw is None:
        return None
    cifmol = cifmol_nw
    pdbid = cifmol.id[0]
    alt_id = cifmol.alt_id
    seq_id_list = []
    seq_cluster_list = []
    for full_chain_id in cifmol.chains.chain_id.value:
        chain_id = full_chain_id.split("_")[0]
        cif_key = f"{pdbid}_{chain_id}_{alt_id}"
        if cif_key not in seq_metadata_map:
            msg = f"Sequence metadata not found for CIF key {cif_key}."
            raise KeyError(msg)
        seqid, seqcluster = seq_metadata_map[cif_key]
        seq_id_list.append(seqid)
        seq_cluster_list.append(seqcluster)
    seq_id_list = np.array(seq_id_list)
    seq_cluster_list = np.array(seq_cluster_list)
    seq_id_list = NodeFeature(seq_id_list)
    seq_cluster_list = NodeFeature(seq_cluster_list)
    cifmol_dict = cifmol.to_dict()
    cifmol_dict["chains"]["nodes"]["seq_id"] = {
        "value": np.array(seq_id_list, dtype=str),
    }
    cifmol_dict["chains"]["nodes"]["cluster_id"] = {
        "value": np.array(seq_cluster_list, dtype=str),
    }
    return cast("dict", CIFMolAttached.from_dict(cifmol_dict).to_dict())


def extract_metadata(
    cifmol_dict: dict[str, dict[str, CIFMol]],
) -> dict[str, dict[str, str]]:
    """Extract resolution and date etc from a CIFMol object."""
    metadata_dict = {}
    NA_types = {
        "polydeoxyribonucleotide",
        "polyribonucleotide",
        "polydeoxyribonucleotide/polyribonucleotide hybrid",
    }
    D_types = {"polypeptide(D)"}

    for cif_key, cifmol_wrapper in cifmol_dict.items():
        cifmol = cifmol_wrapper["cifmol"]
        metadata_dict[cif_key] = {}
        resolution = cifmol.metadata.get("resolution", "NA")
        deposition_date = cifmol.metadata.get("deposition_date", "NA")
        release_date = cifmol.metadata.get("release_date", "NA")
        chain_num = len(cifmol.chains)
        residue_num = len(cifmol.residues)
        atom_num = len(cifmol.atoms)
        including_NA = set(cifmol.chains.entity_type.value).intersection(NA_types)
        including_NA = "Yes" if including_NA else "No"
        including_Dform = set(cifmol.chains.entity_type.value).intersection(D_types)
        including_Dform = "Yes" if including_Dform else "No"
        metadata_dict[cif_key]["resolution"] = str(resolution)
        metadata_dict[cif_key]["deposition_date"] = str(deposition_date)
        metadata_dict[cif_key]["release_date"] = str(release_date)
        metadata_dict[cif_key]["chain_num"] = str(chain_num)
        metadata_dict[cif_key]["residue_num"] = str(residue_num)
        metadata_dict[cif_key]["atom_num"] = str(atom_num)
        metadata_dict[cif_key]["including_NA"] = including_NA
        metadata_dict[cif_key]["including_Dform"] = including_Dform

    return metadata_dict


def classify_seq_clusters(
    raw_fasta_dict: dict[str, str],
    fasta_dict: dict[str, str],
    seqid_map: dict[str, dict[str, str]],
    seqclusters2seqids: dict[str, list[str]],
    db_code: str,
) -> set[str]:
    """Classify sequence clusters based on the provided fasta dictionary and sequence ID map."""
    seqid2seqcluster: dict[str, str] = reverse_dict(seqclusters2seqids)
    classified_clusters: set[str] = set()
    for header in fasta_dict:
        mol_type = header.split("|")[1].strip()
        mol_identifier = mol_type_map.get(mol_type, "X")
        raw_sequence = raw_fasta_dict[header]
        seqid = seqid_map[mol_identifier].get(raw_sequence)
        if seqid is None:
            msg = f"Sequence ID not found for molecule type {mol_identifier} and sequence {raw_sequence}."
            raise KeyError(msg)
        # Absent from this DB's clustering -> its own singleton (see build_seq_metadata_map).
        seqcluster = seqid2seqcluster.get(seqid, cluster_id(db_code, seqid))
        classified_clusters.add(seqcluster)
    return classified_clusters


def load_fasta(
    fasta_path: Path,
) -> dict[str, str]:
    """Load a FASTA file and return its contents as a dictionary."""
    fasta_dict = {}
    with fasta_path.open("r", encoding="utf-8") as f:
        lines = f.readlines()
        current_header = None
        for _line in lines:
            line = _line.strip()
            if line.startswith(">"):
                current_header = line[1:]  # Remove the '>' character
                fasta_dict[current_header] = ""
            elif current_header is not None:
                fasta_dict[current_header] += line
            else:
                msg = "FASTA format error: sequence data found before any header."
                raise ValueError(msg)
    return fasta_dict


def chunk_protein_seqs(
    seqid2seq: dict[str, list[str]],
    chunk_size: int = 5000,
) -> list[dict]:
    """Group protein sequences into chunks that share one MMseqs2 search.

    MMseqs2 pays its cost per search (loading the target DB), not per query, so a
    per-sequence work item would be slower than HHblits rather than faster. Chunking is
    what makes the batched search worth using -- see transforms/msa_mmseqs.py.
    """
    flat = extract_protein_seqs(seqid2seq)
    chunks: list[dict] = []
    for start in range(0, len(flat), chunk_size):
        piece = flat[start : start + chunk_size]
        chunks.append(
            {
                "chunk": [(item["seqid"], item["sequence"]) for item in piece],
                "chunk_index": len(chunks),
            },
        )
    return chunks


def extract_protein_seqs(
    seqid2seq: dict[str, list[str]],
    *,
    remove_unknown: bool = True,
) -> list[dict]:
    """Extract protein sequences from the seqid2seq mapping."""
    protein_seqs = []
    for seqid, seqs in seqid2seq.items():
        if len(seqs) == 0:
            msg = f"No sequence found for sequence ID {seqid}."
            raise ValueError(msg)
        if len(seqs) > 1:
            msg = f"Multiple sequences found for sequence ID {seqid}."
            raise ValueError(msg)
        seq = seqs[0]
        mol_identifier = seqid[0]
        if mol_identifier == "P":
            if remove_unknown:
                is_unknown = all(aa == "X" for aa in seq)
                if is_unknown:
                    continue
            protein_seqs.append({"seqid": seqid, "sequence": seq})
    return protein_seqs


def extract_rna_seqs(
    seqid2seq: dict[str, list[str]],
    *,
    remove_unknown: bool = True,
) -> list[dict]:
    """Extract RNA sequences (``R``-prefixed seq ids) from the seqid2seq mapping.

    Mirror of :func:`extract_protein_seqs` for the RNA moltype so the same
    ``parallel-run`` split/search plumbing can drive RNA MSA generation.
    """
    rna_seqs = []
    for seqid, seqs in seqid2seq.items():
        if len(seqs) == 0:
            msg = f"No sequence found for sequence ID {seqid}."
            raise ValueError(msg)
        if len(seqs) > 1:
            msg = f"Multiple sequences found for sequence ID {seqid}."
            raise ValueError(msg)
        seq = seqs[0]
        mol_identifier = seqid[0]
        if mol_identifier == "R":
            if remove_unknown and all(nt == "N" for nt in seq):
                continue
            rna_seqs.append({"seqid": seqid, "sequence": seq})
    return rna_seqs


def _read_pdb_rna_sequences(pdb_fasta_path: Path) -> set[str]:
    """Collect unique polyribonucleotide sequences from a BioMolDB chain FASTA.

    Headers look like ``>100D_A_. | polyribonucleotide | Auth:A``; only chains
    whose moltype field is exactly ``polyribonucleotide`` (pure RNA, not
    DNA/RNA hybrids) are kept.
    """
    rna_seqs: set[str] = set()
    is_rna = False
    chunks: list[str] = []
    with pdb_fasta_path.open("r", encoding="utf-8") as f:
        for line in f:
            if line.startswith(">"):
                if is_rna and chunks:
                    rna_seqs.add("".join(chunks))
                chunks = []
                fields = line.split(" | ")
                is_rna = len(fields) > 1 and fields[1].strip() == "polyribonucleotide"
            elif is_rna:
                chunks.append(line.strip())
        if is_rna and chunks:
            rna_seqs.add("".join(chunks))
    return rna_seqs


def extract_pdb_rna_seqs(
    seqid2seq: dict[str, list[str]],
    *,
    pdb_fasta_path: Path,
    remove_unknown: bool = True,
) -> list[dict]:
    """Extract PDB-only RNA work items (one per unique polyribonucleotide seq).

    ``seqid2seq`` (from the global seq_id_map) also contains RNA-distillation
    sequences, so restrict to sequences that actually appear as
    ``polyribonucleotide`` chains in the PDB chain FASTA, mapping each back to
    its ``R``-prefixed seq id.
    """
    pdb_rna = _read_pdb_rna_sequences(Path(pdb_fasta_path))
    seq_to_id: dict[str, str] = {}
    for seqid, seqs in seqid2seq.items():
        if seqid[:1] == "R" and seqs:
            seq_to_id.setdefault(seqs[0], seqid)
    out: list[dict] = []
    seen: set[str] = set()
    for seq in pdb_rna:
        seqid = seq_to_id.get(seq)
        if seqid is None or seqid in seen:
            continue
        if remove_unknown and all(nt == "N" for nt in seq):
            continue
        seen.add(seqid)
        out.append({"seqid": seqid, "sequence": seq})
    # Shortest first: nhmmer cost scales with query length, so processing short
    # RNAs first clears the bulk of the set quickly and lets every worker focus
    # on the expensive long rRNA tail instead of starving behind it.
    out.sort(key=lambda item: len(item["sequence"]))
    return out


def parse_metadata(
    metadata_path: Path,
) -> dict[str, time.struct_time]:
    """Map PDB id -> the date used for the AF3/OpenFold3 template time cutoff.

    Header-aware: prefers the ``release_date`` column (initial PDB release),
    falling back to ``deposition_date`` for backward compatibility. The id column
    may be a bare ``pdbid`` or a ``cif_id`` (``{pdbid}_...``); either resolves to
    the lowercase PDB id.
    """
    metadata_dict = {}
    with metadata_path.open("r", encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
        col = {name: i for i, name in enumerate(header)}
        id_idx = col.get("pdbid", col.get("cif_id", 0))
        if "release_date" in col:
            date_idx = col["release_date"]
        elif "deposition_date" in col:
            date_idx = col["deposition_date"]
        else:  # legacy positional layout: cif_id, resolution, deposition_date, ...
            date_idx = 2
        for line in f:
            parts = line.rstrip("\n").split("\t")
            pdb_id = parts[id_idx].split("_")[0].lower()
            date_str = parts[date_idx]
            if not date_str:
                continue
            metadata_dict[pdb_id] = time.strptime(date_str, "%Y-%m-%d")
    return metadata_dict


def build_chain2seqid_map(
    seq_metadata_map: dict[str, tuple[str, str]],
) -> dict[str, str]:
    """Map ``{pdbid}_{chain}`` -> seq_id for chain-keyed template builds."""
    chain2seqid: dict[str, str] = {}
    for cif_id, (seqid, _) in seq_metadata_map.items():
        parts = cif_id.split("_")
        chain2seqid[f"{parts[0]}_{parts[1]}"] = seqid
    return chain2seqid


def build_template_metadata_map(
    pdb_id2deposition_date: dict[str, time.struct_time],
    seq_metadata_map: dict[str, tuple[str, str]],
    seqid2seq: dict[str, list[str]],
    signalp_dict: dict[str, tuple[int, int]],
) -> tuple[dict[str, dict[str, str]], dict[str, time.struct_time], dict[str, str]]:
    """Build a metadata map from sequence cluster ID to template metadata."""
    template_metadata_map = {}
    seqid2earliest_date: dict[str, time.struct_time] = {}
    filtered_seqid2seq: dict[str, str] = {}  # remove signal peptide
    for cif_id, (seqid, _) in seq_metadata_map.items():
        pdb_id = cif_id.split("_")[0]
        chain_id = cif_id.split("_")[1]
        deposition_date = pdb_id2deposition_date.get(pdb_id.lower())
        if deposition_date is None:
            msg = f"Deposition date not found for PDB ID {pdb_id}."
            raise KeyError(msg)
        seq = seqid2seq[seqid][0]
        if seqid in signalp_dict:
            seq = seq[signalp_dict[seqid][1] + 1 :]
        template_metadata_map[f"{pdb_id}_{chain_id}"] = {
            "deposition_date": deposition_date,
            "sequence": seq,
        }
        filtered_seqid2seq[seqid] = seq
        if (
            seqid not in seqid2earliest_date
            or deposition_date < seqid2earliest_date[seqid]
        ):
            seqid2earliest_date[seqid] = deposition_date
    return template_metadata_map, seqid2earliest_date, filtered_seqid2seq


def load_chain_templates(path: Path) -> dict[str, list[str]]:
    r"""Chain -> [template ids] from Phase 2's chain_to_templates.tsv.

    Lines for chains with no templates ('chain\t') are skipped, so a chain is
    present only if it has >=1 template.
    """
    out: dict[str, list[str]] = {}
    with Path(path).open() as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2 or not parts[1]:
                continue
            out[parts[0]] = parts[1].split(",")
    return out


def invert_seqid_to_chains(path: Path) -> dict[str, str]:
    """Chain -> seq_id, inverted from Phase 1's seqid_to_chains.tsv."""
    out: dict[str, str] = {}
    with Path(path).open() as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2 or not parts[1]:
                continue
            for chain in parts[1].split(","):
                out[chain] = parts[0]
    return out


def load_seq_tsv(path: Path) -> dict[str, str]:
    r"""Load a two-column '<key>\t<sequence>' TSV into a dict (single value)."""
    out: dict[str, str] = {}
    with Path(path).open() as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2 and parts[1]:
                out[parts[0]] = parts[1]
    return out
