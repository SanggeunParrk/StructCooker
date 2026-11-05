import lmdb
import pytest
import yaml

from structcooker.instructions.interface.config import DEFAULT, load_config
from structcooker.instructions.interface.db import InterfaceDB
from structcooker.instructions.transforms.codecs import to_bytes


def test_default_config_file_exists_and_loads():
    assert DEFAULT.is_file()
    cfg = load_config()
    assert cfg.contact.cutoff == 5.0
    assert cfg.unit.max_polymer_chains == 6
    assert ("protein", "protein") in cfg.unit.allowed_pairs
    assert (
        cfg.unit.max_polymer_chains_in_assembly
        == cfg.raw["assembly"]["max_polymer_chains"]
    )
    assert cfg.unit.weight_bucket == cfg.raw["dedup"]["weight_bucket"]
    assert cfg.min_heavy == 6
    assert (cfg.water_near, cfg.water_bridge) == (4.0, 3.5)
    assert cfg.duplex["min_bp"] == 4
    assert isinstance(cfg.excluded, frozenset)
    assert len(cfg.excluded) > 0


def test_dotted_overrides_replace_nested_values():
    cfg = load_config(overrides={"unit.max_ligands": 9, "contact.cutoff": 6.5})
    assert cfg.unit.max_ligands == 9
    assert cfg.contact.cutoff == 6.5
    assert cfg.raw["unit"]["max_ligands"] == 9


def test_override_of_a_missing_section_raises():
    with pytest.raises(KeyError):
        load_config(overrides={"nosuch.key": 1})


def test_max_forced_polymer_defaults_to_10_when_absent(tmp_path):
    raw = yaml.safe_load(DEFAULT.read_text())
    del raw["unit"]["max_forced_polymer"]
    p = tmp_path / "i.yaml"
    p.write_text(yaml.safe_dump(raw))
    assert load_config(p).unit.max_forced_polymer == 10


def _write_shard(path, records):
    path.parent.mkdir(parents=True, exist_ok=True)
    env = lmdb.open(str(path), map_size=1 << 24)
    with env.begin(write=True) as txn:
        for k, v in records.items():
            txn.put(k.encode(), to_bytes(v))
    env.close()


REC = {"meta": {"unit_id": "1abc_1_1_0", "n": 3}}


def test_get_searches_all_shards_and_honours_the_shard_hint(tmp_path):
    shards = tmp_path / "shards"
    _write_shard(shards / "shard_0.lmdb", {"1abc_1_1_0": REC})
    _write_shard(shards / "shard_1.lmdb", {"2xyz_1_1_0": {"meta": {"n": 7}}})
    db = InterfaceDB(tmp_path)
    assert sorted(db.shards) == [0, 1]
    assert db.get("1abc_1_1_0")["meta"]["n"] == 3
    assert db.get("2xyz_1_1_0", shard=1)["meta"]["n"] == 7
    with pytest.raises(KeyError):
        db.get("2xyz_1_1_0", shard=0)
    with pytest.raises(KeyError):
        db.get("nope")


def test_published_collection_layout_is_preferred_over_plain_shards(tmp_path):
    _write_shard(tmp_path / "shards" / "shard_0.lmdb", {"a": {"meta": {"n": 1}}})
    _write_shard(
        tmp_path / "interface_pdb.lmdb" / "shards" / "shard_12.lmdb",
        {"a": {"meta": {"n": 2}}},
    )
    db = InterfaceDB(tmp_path)
    assert list(db.shards) == [12]
    assert db.get("a")["meta"]["n"] == 2
