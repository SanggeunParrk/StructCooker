from pathlib import Path

from datacooker import RecipeBook

from structcooker.instructions.transforms.template import list_afm_msa_wo_lower

"""List the AFDB entity MSAs to strip into <seq_id>.a3m (one per seq_id, unfinished only)."""

recipe = RecipeBook()

recipe.step(
    outputs=(("data_list", list),),
    instruction=list_afm_msa_wo_lower,
    kwargs={
        "msa_dir": ("msa_dir", Path),
        "afm_msa_seqid_path": ("afm_msa_seqid_path", Path),
        "output_dir": ("output_dir", Path),
    },
)


RECIPE = recipe
TARGETS = ["data_list"]
