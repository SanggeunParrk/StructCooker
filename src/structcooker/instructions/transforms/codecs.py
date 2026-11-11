import json
import uuid
from io import BytesIO
from typing import cast

import numpy as np
from biomol.core.container import FeatureContainer
from biomol.core.index import IndexTable
from zstandard import ZstdCompressor, ZstdDecompressor, frame_content_size

# First bytes of a Re-Glyco delta value (``transforms.delta_codec``). It is not a zstd frame, so a plain
# ``from_bytes`` can tell it from a BioMol value and refuse it by name instead of failing inside zstd.
DELTA_MAGIC = b"BMDELTA1"


def flatten_data(data: dict) -> tuple[dict, dict]:
    """Flatten a nested dictionary."""
    template = {}
    flatten = {}
    for key, value in data.items():
        if isinstance(value, np.ndarray):
            _key = str(uuid.uuid4())
            template[key] = _key
            buffer = BytesIO()
            np.save(buffer, np.ascontiguousarray(value), allow_pickle=True)
            flatten[_key] = buffer.getvalue()

        elif isinstance(value, dict):
            _template, _flatten = flatten_data(value)
            template[key] = _template
            flatten.update(_flatten)

        elif isinstance(value, FeatureContainer):
            _template, _flatten = flatten_data(cast("dict", value.to_dict()))
            template[key] = _template
            flatten.update(_flatten)

        elif isinstance(value, IndexTable):
            _template, _flatten = flatten_data(indextable_to_dict(value))
            template[key] = _template
            flatten.update(_flatten)

        else:
            template[key] = value
    return template, flatten


def to_raw(data_dict: dict) -> bytes:
    """Serialize the container to its uncompressed bytes (what ``to_bytes`` compresses)."""
    template, flattened_data = flatten_data(data_dict)
    header = {
        "template": template,
        "arrays": {key: len(value) for key, value in flattened_data.items()},
    }
    header_bytes = json.dumps(header).encode("utf-8")
    payload = b"".join(flattened_data[key] for key in flattened_data)
    return len(header_bytes).to_bytes(8, "little") + header_bytes + payload


def to_bytes(data_dict: dict, level: int = 6) -> bytes:
    """Serialize the container to zstd-compressed bytes.

    Parameters
    ----------
    level: int, optional
        The compression level for zstd (default is 6).
    """
    return ZstdCompressor(level=level).compress(to_raw(data_dict))


def keep_bytes(value: bytes) -> bytes:
    """Writer serializer for a value that is already stored bytes (a Re-Glyco delta record)."""
    if not isinstance(value, (bytes, bytearray)):
        msg = f"keep_bytes expects bytes, got {type(value).__name__}"
        raise TypeError(msg)
    return bytes(value)


def reconstruct_data(template: dict, flatten: dict) -> dict:
    """Reconstruct a nested dictionary from its flattened form."""
    data = {}
    for key, value in template.items():
        if isinstance(value, str) and value in flatten:
            buffer = BytesIO(flatten[value])
            buffer.seek(0)
            arr = np.load(buffer, allow_pickle=False)
            data[key] = arr
        elif isinstance(value, dict):
            data[key] = reconstruct_data(value, flatten)
        else:
            data[key] = value
    return data


def from_bytes(byte_data: bytes) -> dict:
    """Deserialize the container from zstd-compressed bytes."""
    if bytes(byte_data[: len(DELTA_MAGIC)]) == DELTA_MAGIC:
        msg = (
            "This value is a Re-Glyco delta record: it can only be decoded against its base record. "
            "Read it with structcooker.instructions.transforms.delta_codec.read_record(env, key)."
        )
        raise ValueError(msg)
    if frame_content_size(byte_data) == -1:
        decoder = ZstdDecompressor().decompressobj()
        raw = decoder.decompress(byte_data)
        if not decoder.eof:
            msg = "Truncated Zstandard frame."
            raise ValueError(msg)
    else:
        raw = ZstdDecompressor().decompress(byte_data)
    return from_raw(raw)


def from_raw(raw: bytes) -> dict:
    """Deserialize the container from its uncompressed bytes (what ``to_raw`` writes)."""
    if len(raw) < 8:
        msg = "Truncated BioMol header length."
        raise ValueError(msg)
    hlen = int.from_bytes(raw[:8], "little")
    if hlen > len(raw) - 8:
        msg = "Truncated BioMol header."
        raise ValueError(msg)
    header = json.loads(raw[8 : 8 + hlen].decode("utf-8"))
    payload = raw[8 + hlen :]
    if sum(header["arrays"].values()) != len(payload):
        msg = "BioMol array payload length does not match its header."
        raise ValueError(msg)

    offset = 0
    flatten_data = {}
    for key, ln in header["arrays"].items():
        chunk = payload[offset : offset + ln]
        offset += ln
        flatten_data[key] = chunk

    template_dict = header["template"]
    return reconstruct_data(template_dict, flatten_data)


def read_header(byte_data: bytes) -> dict:
    """Read scalar metadata without allocating or decoding the array payload."""
    with ZstdDecompressor().stream_reader(BytesIO(byte_data)) as source:
        size = source.read(8)
        if len(size) != 8:
            msg = "Truncated BioMol header length."
            raise ValueError(msg)
        length = int.from_bytes(size, "little")
        header = source.read(length)
        if len(header) != length:
            msg = "Truncated BioMol header."
            raise ValueError(msg)
    return json.loads(header)


def indextable_to_dict(index_table: IndexTable) -> dict:
    """Convert an IndexTable instance to a regular dictionary."""
    return {
        "atom_to_res": index_table.atom_to_res.tolist(),
        "res_to_chain": index_table.res_to_chain.tolist(),
        "res_atom_indptr": index_table.res_atom_indptr.tolist(),
        "res_atom_indices": index_table.res_atom_indices.tolist(),
        "chain_res_indptr": index_table.chain_res_indptr.tolist(),
        "chain_res_indices": index_table.chain_res_indices.tolist(),
    }


def to_dict(data: dict) -> dict:
    """Convert all FeatureContainer and IndexTable instances to regular dictionaries."""
    result = {}
    for key, value in data.items():
        if isinstance(value, FeatureContainer):
            result[key] = to_dict(cast("dict", value.to_dict()))
        elif isinstance(value, IndexTable):
            result[key] = indextable_to_dict(value)
        elif isinstance(value, dict):
            result[key] = to_dict(value)
        else:
            result[key] = value
    return result
