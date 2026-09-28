from pathlib import Path

from datacooker import RecipeBook

from structcooker.instructions.transforms.msa_mmseqs import run_mmseqs_msa_search

"""Build a chunked MMseqs2 MSA search Cooker.

One work item = one chunk of queries sharing a single search against uniref30_2302.
Contrast ``msa_search_hhblits``, where one item is one sequence because that is how
HHblits works. Output is ``<seq_id>.a3m`` in the sharded layout db/msa/a3m.yaml ingests,
so the two paths are interchangeable downstream.
"""

recipe = RecipeBook()

recipe.step(
    outputs=(("result", dict),),
    instruction=run_mmseqs_msa_search,
    kwargs={
        "chunk": ("chunk", list),
        "chunk_index": ("chunk_index", int),
        "work_dir": ("work_dir", Path),
        "output_dir": ("output_dir", Path),
        "db_uniref": ("db_uniref", str | Path),
        "mmseqs_bin": ("mmseqs_bin", str | Path),
        "threads": ("threads", int),
    },
)

RECIPE = recipe
TARGETS = ["result"]
