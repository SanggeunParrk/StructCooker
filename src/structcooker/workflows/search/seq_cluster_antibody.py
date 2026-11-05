from pathlib import Path

from datacooker import RecipeBook

from structcooker.instructions.transforms.sequence import (
    antibody_cluster,
    merge_cluster,
    protein_cluster,
    separate_antibody_sequences,
)

"""``seq_cluster`` for a DB of antibodies only (ASD's ASB).

The antibody/protein split reads the role from the fasta headers (``Auth:H`` / ``Auth:L``) instead of SAbDab,
and the CDR-H3/L3 extraction runs ANARCI in batches. Everything after the split is ``seq_cluster``'s.
"""

seq_cluster_antibody_recipe = RecipeBook()

seq_cluster_antibody_recipe.step(
    outputs=(("fasta_path_dict", dict),),
    instruction=separate_antibody_sequences,
    kwargs={
        "tmp_dir": ("tmp_dir", Path),
        "seq_id_map_path": ("seq_id_map_path", Path),
        "fasta_path": ("fasta_path", Path),
    },
)


seq_cluster_antibody_recipe.step(
    outputs=(
        ("antibody_cluster_dict", dict),
        ("failed_ids", set),
    ),
    instruction=antibody_cluster,
    kwargs={
        "tmp_dir": ("tmp_dir", Path),
        "fasta_path_dict": ("fasta_path_dict", Path),
        "anarci_batch_size": ("anarci_batch_size", int),
        "anarci_ncpu": ("anarci_ncpu", int),
        "h3_identity": ("h3_identity", float),
        "l3_identity": ("l3_identity", float),
    },
)

seq_cluster_antibody_recipe.step(
    outputs=(
        ("protein_cluster_dict", dict),
        ("protein_d_cluster_dict", dict),
    ),
    instruction=protein_cluster,
    kwargs={
        "tmp_dir": ("tmp_dir", Path),
        "fasta_path_dict": ("fasta_path_dict", Path),
        "failed_ids": ("failed_ids", set),
        "seq_id_map_path": ("seq_id_map_path", Path),
        "mmseqs2_seq_id": ("mmseqs2_seq_id", float | str),
        "mmseqs2_cov": ("mmseqs2_cov", float | str),
        "mmseqs2_covmode": ("mmseqs2_covmode", str),
        "mmseqs2_clustermode": ("mmseqs2_clustermode", str),
    },
)

seq_cluster_antibody_recipe.step(
    outputs=(("db_code", str),),
    instruction=lambda db_code: db_code,
    kwargs={"db_code": ("db_code", str)},
)

seq_cluster_antibody_recipe.step(
    outputs=(("cluster_dict", dict),),
    instruction=merge_cluster,
    kwargs={
        "fasta_path_dict": ("fasta_path_dict", dict),
        "protein_cluster_dict": ("protein_cluster_dict", dict),
        "protein_d_cluster_dict": ("protein_d_cluster_dict", dict),
        "antibody_cluster_dict": ("antibody_cluster_dict", dict),
    },
)

RECIPE = seq_cluster_antibody_recipe
TARGETS = ["cluster_dict", "db_code"]
