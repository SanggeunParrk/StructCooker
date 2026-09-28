import numpy as np
import pytest
from zstandard import ZstdCompressor, ZstdDecompressor

from structcooker.instructions.transforms.codecs import from_bytes, to_bytes


def test_legacy_frame_without_content_size():
    source = {"msa_dict": {"array": np.arange(20).reshape(4, 5)}}
    normal = to_bytes(source)
    raw = ZstdDecompressor().decompress(normal)
    streaming = ZstdCompressor(write_content_size=False).compress(raw)
    np.testing.assert_array_equal(from_bytes(streaming)["msa_dict"]["array"], source["msa_dict"]["array"])
    with pytest.raises(ValueError, match="Truncated"):
        from_bytes(streaming[:-3])


def test_incomplete_array_payload_is_rejected():
    raw = ZstdDecompressor().decompress(to_bytes({"array": np.arange(10)}))
    with pytest.raises(ValueError, match="payload length"):
        from_bytes(ZstdCompressor().compress(raw[:-1]))
