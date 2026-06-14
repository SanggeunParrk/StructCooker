from pathlib import Path

from datacooker import RecipeBook

from structcooker.instructions.transforms.template import run_hmmbuild, run_hmmsearch

"""Build a FASTA->MSA (SignalP + HHblits/HHfilter) Cooker."""

recipe = RecipeBook()

recipe.step(
    outputs=(("hmmbuild_results", dict),),
    instruction=run_hmmbuild,
    kwargs={
        "input_a3m_path": ("input_a3m_path", Path),
        "hmm_path": ("output_path", Path),
    },
)


recipe.step(
    outputs=(("hmmsearch_results", dict),),
    instruction=run_hmmsearch,
    kwargs={
        "output_dir": ("hmm_output_dir", Path),
        "hmm_path": ("output_path", Path),
        "fasta_path": ("fasta_path", Path),
    },
)


RECIPE = recipe
TARGETS = ["hmmsearch_results", "hmmbuild_results"]
