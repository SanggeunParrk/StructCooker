import logging
import os
import tempfile
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
