from pathlib import Path

from datacooker import RecipeBook

from structcooker.instructions.readers.sequence import load_fasta
from structcooker.instructions.transforms.metadata import (
    build_seq_metadata_map,
    load_tsv,
    seq_metadata_map_to_lmdb,
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
    outputs=(("seq_metadata_map_dict", dict),),
    instruction=build_seq_metadata_map,
    kwargs={
        "raw_fasta_dict": ("raw_fasta_dict", dict),
        "seqid2seq": ("seqid2seq", dict),
        "seqclusters2seqids": ("seqclusters2seqids", dict),
        "db_code": ("db_code", str),
    },
)

# Wrap the multi-GB dict as an on-disk LMDB so Ray workers mmap one shared copy per
# node instead of each deserializing their own (the per-worker metadata OOM). The dict
# is built once on the shard driver, written here, and only the LMDB path travels to
# workers. attach_metadata uses it exactly like the dict (``in`` / ``[]``).
recipe.step(
    outputs=(("seq_metadata_map", object),),
    instruction=seq_metadata_map_to_lmdb,
    kwargs={
        "seq_metadata_map": ("seq_metadata_map_dict", dict),
    },
)

# Only seq_metadata_map is shipped to the rebuild workers (attach_metadata consumes
# just that -- a small cif_id -> (seq_id, cluster) map over PDB chains). The heavy
# intermediates seqid2seq (the 6.6 GB seq_id_map) and seqclusters2seqids are computed
# on the driver to BUILD it, but must NOT be targets: as targets they rode into every
# Ray worker's metadata copy (6.6 GB x 112 workers -> hundreds of GB -> OOM, even on
# "tiny" items, since the footprint was per-worker metadata, not per-item).
RECIPE = recipe
TARGETS = ["seq_metadata_map"]
