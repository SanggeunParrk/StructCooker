import numpy as np

from structcooker.instructions.transforms.codecs import read_header, to_bytes
from structcooker.instructions.transforms.template import (
    adapt_cif_record_to_chain_inputs,
)


def _assembly():
    return {"atoms": {"nodes": {"id": {"value": np.array(["CA"])}, "occupancy": {"value": np.array([1.])}}},
            "chains": {"nodes": {"chain_id": {"value": np.array(["A_1"])}}}}


def test_equal_score_model_choice_is_numeric_and_order_independent():
    for order in [("1_10_.", "1_2_."), ("1_2_.", "1_10_.")]:
        inputs=adapt_cif_record_to_chain_inputs({"assembly_dict":{key:_assembly() for key in order}})
        assert inputs["A"]["cif_key"] == "1_2_."


def test_header_reads_metadata_without_loading_arrays(monkeypatch):
    payload=to_bytes({"metadata":{"model_id":"10"},"array":np.arange(100)})
    def fail(*_args, **_kwargs):
        message = "array decoder must not run"
        raise AssertionError(message)
    monkeypatch.setattr(np,"load",fail)
    assert read_header(payload)["template"]["metadata"] == {"model_id":"10"}


def test_chain_model_cache_is_shared_only_within_one_record(monkeypatch):
    from types import SimpleNamespace

    from structcooker.instructions.transforms.template import extract_selected_chain

    built = []
    class Selection:
        chain_id = SimpleNamespace(value=np.array(["A_1"]))
        def __getitem__(self, _mask):
            return self
        def extract(self):
            return SimpleNamespace(to_dict=lambda: {"ok": True})
    def construct(data):
        built.append(data)
        return SimpleNamespace(chains=Selection())
    monkeypatch.setattr("structcooker.instructions.transforms.template.CIFMol.from_dict", construct)
    record = {"assembly_dict": {"1_1_.": _assembly(), "1_2_.": _assembly()}, "metadata_dict": {"id": ["test"]}}
    prepared = adapt_cif_record_to_chain_inputs(record)
    cache_state = prepared["A"]["model_cache"]
    for _ in range(3):
        assert extract_selected_chain(record, "A", "1_1_.", model_cache=cache_state) == {"ok": True}
    assert len(built) == 1
    extract_selected_chain(record, "A", "1_2_.", model_cache=cache_state)
    assert len(built) == 2
    fresh = {**record}
    extract_selected_chain(fresh, "A", "1_2_.", model_cache=cache_state)
    assert len(built) == 3
    assert adapt_cif_record_to_chain_inputs(fresh)["A"]["model_cache"] is not cache_state


def test_atom_cap_drops_oversized_assemblies_not_the_entry(monkeypatch):
    import structcooker.instructions.transforms.template as tmpl

    def assembly(n_atoms, occupancy):
        return {"atoms": {"nodes": {"id": {"value": np.array(["CA"] * n_atoms)},
                                    "occupancy": {"value": np.full(n_atoms, occupancy)}}},
                "chains": {"nodes": {"chain_id": {"value": np.array(["A_1"])}}}}

    monkeypatch.setattr(tmpl, "CIFCHAIN_MAX_ATOMS", 10)
    # A capsid-like entry: the biological assembly is over the cap, the asymmetric unit is not.
    # The chain must come from the small assembly even though the big one scores higher.
    record = {"assembly_dict": {"1_1_.": assembly(50, 1.0), "2_1_.": assembly(5, 1.0)}}
    assert tmpl.adapt_cif_record_to_chain_inputs(record)["A"]["cif_key"] == "2_1_."
    # Only when every assembly is over the cap is the entry skipped.
    assert tmpl.adapt_cif_record_to_chain_inputs({"assembly_dict": {"1_1_.": assembly(50, 1.0)}}) == {}
