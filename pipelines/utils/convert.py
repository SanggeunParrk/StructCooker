import json
import uuid
from io import BytesIO
from typing import cast

import numpy as np
from biomol.core.container import FeatureContainer
from biomol.core.index import IndexTable
from zstandard import ZstdCompressor, ZstdDecompressor


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


def to_bytes(data_dict: dict, level: int = 6) -> bytes:
    """Serialize the container to zstd-compressed bytes.

    Parameters
    ----------
    level: int, optional
        The compression level for zstd (default is 6).
    """
    template, flattened_data = flatten_data(data_dict)
    header = {
        "template": template,
        "arrays": {key: len(value) for key, value in flattened_data.items()},
    }
    header_bytes = json.dumps(header).encode("utf-8")
    payload = b"".join(flattened_data[key] for key in flattened_data)
    raw = len(header_bytes).to_bytes(8, "little") + header_bytes + payload
    return ZstdCompressor(level=level).compress(raw)


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
    raw = ZstdDecompressor().decompress(byte_data)
    hlen = int.from_bytes(raw[:8], "little")
    header = json.loads(raw[8 : 8 + hlen].decode("utf-8"))
    payload = raw[8 + hlen :]

    offset = 0
    flatten_data = {}
    for key, ln in header["arrays"].items():
        chunk = payload[offset : offset + ln]
        offset += ln
        flatten_data[key] = chunk

    template_dict = header["template"]
    return reconstruct_data(template_dict, flatten_data)


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
