import json

import numpy as np

from structcooker.instructions.readers import openfold


def test_explicit_structure_recovery_keeps_entry_identity(tmp_path, monkeypatch):
    original = tmp_path / "MGYP1" / "structure.npz"
    original.parent.mkdir()
    original.write_bytes(b"broken original")
    recovered = tmp_path / "recovered.npz"
    np.savez(recovered, chain_id=np.array(["A", "A", "A"]),
             res_id=np.array([1, 1, 2]), ins_code=np.array(["", "", ""]))
    template = original.with_name("template.npz")
    np.savez(template)
    manifest = tmp_path / "overrides.json"
    manifest.write_text(json.dumps({str(original): str(recovered)}))
    monkeypatch.setenv("OPENFOLD_STRUCTURE_OVERRIDES", str(manifest))
    monkeypatch.setattr(openfold, "_STRUCTURE_OVERRIDES", None)
    assert openfold.get_openfold_template_data(template)["query_len"] == 2
    assert openfold.get_openfold_structure_data(original)["entry_id"] == "MGYP1"
    assert original.read_bytes() == b"broken original"


def test_run_scoped_sequence_map(tmp_path, monkeypatch):
    mapping = tmp_path / "ids.tsv"
    mapping.write_text("MGYP1\tP123\n")
    monkeypatch.setenv("MONOMER_SEQID_MAP", str(mapping))
    monkeypatch.setattr(openfold, "_MONOMER_SEQID", None)
    assert openfold.openfold_seqid_key(tmp_path / "MGYP1" / "alignment.npz") == "P123"


def test_empty_npz_msa_is_rejected_at_reader(tmp_path):
    import numpy as np
    import pytest

    from structcooker.instructions.readers.openfold import get_openfold_msa_data

    path = tmp_path / "alignment.npz"
    np.savez(path)
    with pytest.raises(ValueError, match="no MSA sources"):
        get_openfold_msa_data(path)
