from pathlib import Path

from datacooker import RecipeBook

from structcooker.instructions.readers.sequence import load_fasta
from structcooker.instructions.transforms.metadata import (
    build_seq_metadata_map,
    build_template_metadata_map,
    load_signalp,
    load_tsv,
    parse_metadata,
)

"""Rebuild a CIF lmdb to train AF3"""

recipe = RecipeBook()

recipe.step(
    outputs=(("pdb_id2deposition_date", dict),),
    instruction=parse_metadata,
    kwargs={
        "metadata_path": ("metadata_path", Path),
    },
)


recipe.step(
    outputs=(("raw_fasta_dict", dict),),
    instruction=load_fasta,
    kwargs={
        "fasta_path": ("raw_fasta_path", str | Path),
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
    outputs=(("seq_metadata_map", dict),),
    instruction=build_seq_metadata_map,
    kwargs={
        "raw_fasta_dict": ("raw_fasta_dict", dict),
        "seqid2seq": ("seqid2seq", dict),
        "seqclusters2seqids": ("seqclusters2seqids", dict),
    },
)


recipe.step(
    outputs=(("signalp_dict", dict),),
    instruction=load_signalp,
    kwargs={
        "signalp_dir": ("signalp_dir", Path | None),
    },
)


recipe.step(
    outputs=(
        ("template_metadata_map", dict),
        ("seqid2earliest_date", dict),
        ("filtered_seqid2seq", dict),
    ),
    instruction=build_template_metadata_map,
    kwargs={
        "pdb_id2deposition_date": ("pdb_id2deposition_date", dict),
        "seq_metadata_map": ("seq_metadata_map", dict),
        "seqid2seq": ("seqid2seq", dict),
        "signalp_dict": ("signalp_dict", dict | None),
    },
)

RECIPE = recipe
TARGETS = ["template_metadata_map", "seqid2earliest_date", "filtered_seqid2seq"]
