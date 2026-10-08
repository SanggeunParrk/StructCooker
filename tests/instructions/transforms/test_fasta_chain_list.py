from structcooker.instructions.transforms.filelists import fasta_chain_list


def test_each_chain_listed_once_across_altlocs(tmp_path):
    fasta = tmp_path / "cif_pdb.fasta"
    fasta.write_text(
        ">1J4R_A_B | polypeptide(L) | Auth:A\nMK\n"
        ">1J4R_D_B | non-polymer | Auth:A\n(HOH)\n"
        ">1J4R_A_A | polypeptide(L) | Auth:A\nMK\n"
        ">1J4R_B_. | polypeptide(L) | Auth:B\nGG\n"
        ">1J4R_A_. | polypeptide(L) | Auth:A\nMK\n",
    )
    out = tmp_path / "chains.txt"
    fasta_chain_list(fasta, out)
    assert out.read_text().split() == ["1J4R_A", "1J4R_B"]
