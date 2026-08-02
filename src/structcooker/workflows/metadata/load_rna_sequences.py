from pathlib import Path

from datacooker import RecipeBook

from structcooker.instructions.transforms.metadata import (
    extract_pdb_rna_seqs,
    load_tsv,
)

"""Split recipe: expand PDB polyribonucleotide chains into per-sequence work items."""

recipe = RecipeBook()

recipe.step(
    outputs=(("seqid2seq", dict),),
    instruction=load_tsv,
    kwargs={
        "tsv_file_path": ("seq_id_map_path", Path),
    },
    params={
        "split_by_comma": False,
    },
)


recipe.step(
    outputs=(("data_list", list),),
    instruction=extract_pdb_rna_seqs,
    kwargs={
        "seqid2seq": ("seqid2seq", dict),
        "pdb_fasta_path": ("pdb_fasta_path", Path),
    },
)


RECIPE = recipe
TARGETS = ["data_list"]
