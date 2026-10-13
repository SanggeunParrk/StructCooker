"""glob_file_list: lists several trees; relative lines resolve against the reader's data_dir."""
from datacooker.cli.lmdb import resolve_listed

from structcooker.instructions.transforms.filelists import glob_file_list


def test_relative_list_resolves_after_the_tree_moves(tmp_path):
    root = tmp_path / "materials"
    for rel in ("raw/msa/000/AF-1-msa_v1.a3m.zst", "raw/msa/001/AF-2-msa_v1.a3m.zst",
                "intermediate/gen/AF-3-msa_v1.a3m.zst", "raw/msa/000/notes.txt"):
        (root / rel).parent.mkdir(parents=True, exist_ok=True)
        (root / rel).write_text(rel)
    sources = [{"dir": str(root / "raw/msa"), "glob": "AF-*", "recursive": True},
               {"dir": str(root / "intermediate/gen"), "glob": "AF-*"}]
    out = tmp_path / "list.txt"
    glob_file_list(sources, out, relative_to=root)
    lines = out.read_text().split()
    assert lines == ["intermediate/gen/AF-3-msa_v1.a3m.zst", "raw/msa/000/AF-1-msa_v1.a3m.zst",
                     "raw/msa/001/AF-2-msa_v1.a3m.zst"]
    moved = root.rename(tmp_path / "renamed")
    assert all(resolve_listed(line, moved).read_text() == line for line in lines)
    glob_file_list([{"dir": str(moved / "raw/msa"), "glob": "AF-*", "recursive": True}], out)   # default: absolute
    assert out.read_text().split()[0] == str(moved / "raw/msa/000/AF-1-msa_v1.a3m.zst")
