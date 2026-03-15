import os
from pathlib import Path
import lmdb

import json
from collections.abc import Mapping
from io import BytesIO
from typing import Any

import numpy as np
from zstandard import ZstdCompressor, ZstdDecompressor


def to_bytes(data: Mapping[str, Any], level: int = 6) -> bytes:
    """Serialize a dictionary containing NumPy arrays to zstd-compressed bytes.

    Parameters
    ----------
    data : Mapping[str, Any]
        A Mapping containing the NumPy arrays and other data to serialize.
    level : int, optional
        The compression level for zstd (default is 6).
    """

    def _flatten_data(
        data: Mapping[str, Any],
    ) -> tuple[dict[str, Any], dict[str, Any]]:
        template = {}
        flatten = {}
        for key, value in data.items():
            if isinstance(value, np.ndarray):
                _key = str(id(value))
                template[key] = _key
                buffer = BytesIO()
                np.save(buffer, np.ascontiguousarray(value), allow_pickle=False)
                flatten[_key] = buffer.getbuffer()
            elif isinstance(value, dict):
                _template, _flatten = _flatten_data(value)
                template[key] = _template
                flatten.update(_flatten)
            else:
                template[key] = value
        return template, flatten

    template, flatten_data = _flatten_data(data)
    header = {
        "template": template,
        "arrays": {key: len(value) for key, value in flatten_data.items()},
    }
    header_bytes = json.dumps(header).encode("utf-8")
    output = BytesIO()
    with ZstdCompressor(level=level).stream_writer(output, closefd=False) as writer:
        writer.write(len(header_bytes).to_bytes(8, "little"))
        writer.write(header_bytes)
        for key in flatten_data:
            writer.write(flatten_data[key])
    return output.getvalue()


def load_bytes(byte_data: bytes) -> Mapping[str, Any]:
    """Deserialize zstd-compressed bytes back into a dictionary."""

    def _reconstruct_data(
        template: dict[str, Any],
        flatten: dict[str, Any],
    ) -> dict[str, Any]:
        data = {}
        for key, value in template.items():
            if isinstance(value, str) and value in flatten:
                buffer = BytesIO(flatten[value])
                arr = np.load(buffer, allow_pickle=False)
                data[key] = arr
            elif isinstance(value, dict):
                data[key] = _reconstruct_data(value, flatten)
            else:
                data[key] = value
        return data

    with ZstdDecompressor().stream_reader(BytesIO(byte_data)) as reader:
        raw = reader.read()
    hlen = int.from_bytes(raw[:8], "little")
    header = json.loads(raw[8 : 8 + hlen].decode("utf-8"))
    payload = memoryview(raw)[8 + hlen :]

    offset = 0
    flatten_data = {}
    for key, ln in header["arrays"].items():
        chunk = payload[offset : offset + ln]
        offset += ln
        flatten_data[key] = chunk

    template_dict = header["template"]
    return _reconstruct_data(template_dict, flatten_data)


def cp_dir(src:Path, dst:Path) -> None:
    if not dst.exists():
        dst.mkdir(parents=True)
    os.system(f"rsync -av --exclude='*.hhr' --exclude='*.atab' {src}/ {dst}/")

def load_seq_id_map(tsv_path:Path) -> dict:
    seq_id_map = {}
    with open(tsv_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue 
            seq_id, seq = line.split("\t")
            if not seq_id.startswith("P"):
                continue
            if seq in seq_id_map:
                assert seq_id_map[seq] == seq_id, f"Duplicate sequence found: {seq}"
                continue
            if set(seq) == {"X"}:
                continue
            if len(seq) < 16:
                continue
            seq_id_map[seq] = seq_id
    return seq_id_map

def load_query_seq_from_a3m(a3m_path:Path) -> str:
    with open(a3m_path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                continue
            return line
    return ""

def main(seq_id_map_path: Path, old_msa_path: Path, new_msa_path: Path) -> None:
    seq_id_map = load_seq_id_map(seq_id_map_path)
    subdir_list = [d for d in old_msa_path.iterdir() if d.is_dir()]
    not_found_count_list = []
    found_list = []
    for subdir in subdir_list:
        old_msa_file = subdir / f"t000_msa0.a3m"
        query_seq = load_query_seq_from_a3m(old_msa_file)
        if query_seq not in seq_id_map:
            # print(f"Query sequence not found in map: {query_seq}")
            not_found_count_list.append(query_seq) # due to signal peptide
            continue
    
        seq_id = seq_id_map[query_seq]
        new_subdir = new_msa_path / seq_id[:4] / seq_id
        new_subdir.mkdir(parents=True, exist_ok=True)
        cp_dir(subdir, new_subdir)
        print(f"Copied {subdir} to {new_subdir}")
        found_list.append(query_seq)
    breakpoint()

def len_lmdb(lmdb_path: Path) -> int:
    env = lmdb.open(
        str(lmdb_path),
        readonly=True,
        lock=False,
        readahead=False,
        max_readers=1
    )

    with env.begin() as txn:
        stat = txn.stat()
        return stat["entries"]
    
def load_item_from_lmdb(lmdb_path:Path, key_idx: int) -> tuple:
    env = lmdb.open(
        str(lmdb_path),
        readonly=True,
        lock=False,
        readahead=False,
        max_readers=1
    )

    with env.begin() as txn:
        cursor = txn.cursor()
        for idx, (key, value) in enumerate(cursor):
            if idx == key_idx:
                value = load_bytes(bytes(value))
                msa_depth = value['msa_container']['residue_container']['nodes']['sequences']['value'].shape[1]
                return key.decode(), value, msa_depth
    return None, None, None

def num_found(msa_db_path: Path) -> int:
    total_count = len_lmdb(msa_db_path)
    found_count = 0
    for idx in range(total_count):
        key, value, msa_depth = load_item_from_lmdb(msa_db_path, idx)
        if key is None:
            continue
        if msa_depth > 0:
            found_count += 1
    return found_count

if __name__ == "__main__":
    seq_id_map_path = Path("/data/psk6950/BioMolDB_20260224/metadata/seq_id_map.tsv")
    old_msa_path = Path("/data/psk6950/MSA/PDB_2024Mar18/new_hash_a3m")
    new_msa_path = Path("/data/psk6950/msa/")	
    main(seq_id_map_path, old_msa_path, new_msa_path)
