"""Clustering a sequence type the build does not contain yields empty clusters, not an error."""
from structcooker.instructions.transforms.sequence import (
    parse_cdhit_cluster,
    parse_mmseqs2_cluster,
    run_CDHIT,
    run_mmseqs2,
)


def test_mmseqs_on_an_empty_fasta_gives_no_clusters(tmp_path):
    fasta = tmp_path / "protein_D.fasta"
    fasta.write_text("")
    run_mmseqs2(fasta, tmp_path / "mmseqs2_protein_D")
    assert parse_mmseqs2_cluster(tmp_path / "mmseqs2_protein_D_cluster.tsv") == {}


def test_cdhit_on_an_empty_fasta_gives_no_clusters(tmp_path):
    fasta = tmp_path / "ab_H3.fasta"
    fasta.write_text("")
    run_CDHIT(fasta, tmp_path / "cdhit_ab_H_clustered.dat")
    assert parse_cdhit_cluster(tmp_path / "cdhit_ab_H_clustered.dat.clstr") == {}
