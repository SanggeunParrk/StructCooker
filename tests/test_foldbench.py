import json

import lmdb
import pytest

from structcooker import foldbench as fb


def make_db(path, items):
    path.mkdir(parents=True)
    with lmdb.open(str(path), map_size=1 << 24) as env, env.begin(write=True) as txn:
        for k, v in items.items():
            txn.put(k.encode(), v)
    return path


def read_all(path):
    with lmdb.open(str(path), readonly=True, lock=False) as env, env.begin() as txn:
        return {bytes(k).decode(): bytes(v) for k, v in txn.cursor()}


def test_lmdb_key_pads_short_ids_and_keeps_full_keys():
    assert fb.lmdb_key("R0021166") == "R00000000000000021166"
    assert fb.lmdb_key("P0300001") == "P00000000000000300001"
    assert fb.lmdb_key("R00000000000000021166") == "R00000000000000021166"
    assert len(fb.lmdb_key("D0016738")) == fb.KEY_WIDTH


def make_package(root):
    meta = root / "meta"
    meta.mkdir(parents=True)
    (meta / "seqids_protein.txt").write_text("P0000001\nP0000002\n")
    (meta / "seqids_rna.txt").write_text("R00000000000000000007\n")
    (meta / "chain_seqid_map.tsv").write_text(
        "target\tchain_id\tentity_type\tseq_id\tn_res\ttrim_offset\tn_res_used\tsequence\n"
        "t1-assembly1\tA\tpolypeptide(L)\tP0000001\t3\t0\t3\t(ACE)AGK\n"
        "t1-assembly1\tB\tpolyribonucleotide\tR0000007\t3\t0\t3\tACG\n"
        "t1-assembly1\tC\tnon-polymer\t\t1\t0\t1\t(ZN)\n"
        "t1-assembly1\tD\tbranched\t\t2\t0\t2\t(NAG)(NAG)|(1,2)\n"
    )
    inputs = root / "inputs"
    (inputs / "t1-assembly1").mkdir(parents=True)
    (inputs / "index.txt").write_text("t1-assembly1\n")
    (inputs / "t1-assembly1" / "data.yaml").write_text(
        "a3m_key:\n  '0': P00000000000000000001\n  '1': R00000000000000000007\n  '2': D00000000000000020007\n"
        "  '3': P00000000000000300001\n"
        "template:\n  '0': t1xx_A\n  '1': t1xx_Z\n"
    )
    return root


def test_load_keys_reads_lists_inline_codes_and_template_keys(tmp_path):
    keys = fb.load_keys(make_package(tmp_path / "pkg"))
    # the data.yaml lookups are authoritative: they name the from-scratch (edge) id and the DNA key that seqids_*.txt do not
    assert keys.protein == ["P00000000000000000001", "P00000000000000300001"]
    assert keys.rna == ["R00000000000000000007"]
    assert keys.other == ["D00000000000000020007"]
    assert keys.ccd == ["ACE", "NAG", "ZN"]  # the glycosidic '|(1,2)' is not a component
    assert keys.template == ["t1xx_A", "t1xx_Z"]


def test_copy_keys_is_byte_exact_reports_missing_and_uses_the_fallback(tmp_path):
    primary = make_db(tmp_path / "primary", {"a": b"1", "b": b"22", "z": b"zz"})
    fallback = make_db(tmp_path / "fallback", {"c": b"333", "a": b"WRONG"})
    rep = fb.copy_keys(tmp_path / "out", ["a", "b", "c", "nope"], [primary, fallback], name="x")
    assert read_all(tmp_path / "out") == {"a": b"1", "b": b"22", "c": b"333"}  # primary wins; 'z' not asked for
    assert (rep.requested, rep.written, rep.from_fallback, rep.missing) == (4, 3, 1, ["nope"])


def test_copy_keys_never_replaces_an_existing_db(tmp_path):
    src = make_db(tmp_path / "src", {"a": b"1"})
    make_db(tmp_path / "out", {"old": b"x"})
    with pytest.raises(FileExistsError):
        fb.copy_keys(tmp_path / "out", ["a"], [src])
    assert read_all(tmp_path / "out") == {"old": b"x"}


def test_build_then_verify_resolves_every_lookup(tmp_path):
    pkg = make_package(tmp_path / "pkg")
    rel = tmp_path / "rel" / "lmdb" / "pdb"
    make_db(rel / "msa" / "a3m.lmdb", {"P00000000000000000001": b"m1", "P00000000000000000002": b"m2"})
    make_db(rel / "msa" / "msa_rna_pdb_shard0.lmdb", {"R00000000000000000007": b"r7"})
    make_db(rel / "template" / "template.lmdb", {"t1xx_A": b"t"})
    edge = make_db(tmp_path / "edge_a3m", {"P00000000000000300001": b"e1"})
    ccd = make_db(tmp_path / "CCD_x.lmdb", {"ZN": b"z", "NAG": b"n"})
    out = tmp_path / "out"
    out.mkdir()
    summary = fb.build(pkg, out, fb.Sources.from_release(tmp_path / "rel", ccd=ccd, fallback=[edge]))
    assert json.loads((out / "foldbench_build.json").read_text())["keys"]["rna"] == 1
    copies = {c["name"]: c for c in summary["copies"]}
    assert copies["protein_msa"]["written"] == 2 and copies["protein_msa"]["from_fallback"] == 1
    assert copies["template"]["missing"] == ["t1xx_Z"] and copies["ccd"]["missing"] == ["ACE"]
    v = fb.verify(pkg, out, decode_sample=0)
    assert v["msa_slots"] == {"found": 3, "expected_missing": 1, "missing": 0}  # the D key has no MSA by design
    assert v["template_slots"] == {"found": 1, "missing": 1}
