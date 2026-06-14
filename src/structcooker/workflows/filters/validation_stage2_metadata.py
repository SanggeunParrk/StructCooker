from pathlib import Path

from datacooker import RecipeBook

from structcooker.instructions.readers.sequence import load_fasta
from structcooker.instructions.transforms.filtering import (
    filter_valid_2_clusters,
)
from structcooker.instructions.transforms.metadata import (
    build_seqid_map,
    classify_seq_clusters,
    load_tsv,
)

"""Rebuild a CIF lmdb to train AF3"""

recipe = RecipeBook()

recipe.step(
    outputs=(("raw_fasta_dict", dict),),
    instruction=load_fasta,
    kwargs={
        "fasta_path": ("raw_fasta_path", str | Path),
    },
)

recipe.step(
    outputs=(("train_fasta_dict", dict),),
    instruction=load_fasta,
    kwargs={
        "fasta_path": ("train_fasta_path", str | Path),
    },
)

recipe.step(
    outputs=(("valid_1_fasta_dict", dict),),
    instruction=load_fasta,
    kwargs={
        "fasta_path": ("valid_1_fasta_path", str | Path),
    },
)


recipe.step(
    outputs=(("seqid2seq", dict),),
    instruction=load_tsv,
    kwargs={
        "tsv_file_path": ("seqid2seq_path", Path),
    },
    params={
        "split_by_comma": False,
    },
)

recipe.step(
    outputs=(("seqclusters2seqids", dict),),
    instruction=load_tsv,
    kwargs={
        "tsv_file_path": ("seqcluster_path", Path),
    },
    params={
        "split_by_comma": True,
    },
)

recipe.step(
    outputs=(("interacting_seq_clusters", dict),),
    instruction=load_tsv,
    kwargs={
        "tsv_file_path": ("interacting_seq_clusters_path", Path),
    },
    params={
        "split_by_comma": False,
    },
)

recipe.step(
    outputs=(("seqid_map", dict),),
    instruction=build_seqid_map,
    kwargs={
        "seqid2seq": ("seqid2seq", dict),
    },
)

recipe.step(
    outputs=(("train_clusters", set),),
    instruction=classify_seq_clusters,
    kwargs={
        "raw_fasta_dict": ("raw_fasta_dict", dict),
        "fasta_dict": ("train_fasta_dict", dict),
        "seqid_map": ("seqid_map", dict),
        "seqclusters2seqids": ("seqclusters2seqids", dict),
    },
)

recipe.step(
    outputs=(("valid_1_clusters", set),),
    instruction=classify_seq_clusters,
    kwargs={
        "raw_fasta_dict": ("raw_fasta_dict", dict),
        "fasta_dict": ("valid_1_fasta_dict", dict),
        "seqid_map": ("seqid_map", dict),
        "seqclusters2seqids": ("seqclusters2seqids", dict),
    },
)

recipe.step(
    outputs=(("filtered_valid_2_clusters", set),),
    instruction=filter_valid_2_clusters,
    kwargs={
        "train_clusters": ("train_clusters", set),
        "valid_1_clusters": ("valid_1_clusters", set),
        "interacting_seq_clusters": ("interacting_seq_clusters", dict),
    },
)

RECIPE = recipe
TARGETS = ["filtered_valid_2_clusters"]
