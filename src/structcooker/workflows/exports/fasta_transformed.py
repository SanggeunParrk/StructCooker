from datacooker import RecipeBook

from structcooker.instructions.transforms.sequence import build_fasta
from structcooker.mols import CIFMol

"""Build fasta from an already-transformed CIF lmdb.

Unlike ``exports/fasta.py`` (which converts raw ``assembly_dict``/``metadata_dict``
via ``convert_to_cifmol_dict``), this variant expects the reader ``adapter`` to
have already produced the ``{cif_key: {"cifmol": CIFMol}}`` structure — i.e. use
``convert_to_cifmol_transformed`` in the yaml. This matches LMDBs written by a
``rebuild`` step (e.g. the 20210930 training filter output).
"""

fasta_recipe = RecipeBook()

fasta_recipe.step(
    outputs=(("fasta", str),),
    instruction=build_fasta,
    kwargs={
        "cifmol_dict": ("db_data", dict[str, dict[str, CIFMol]] | None),
    },
)

RECIPE = fasta_recipe
TARGETS = ["fasta"]
