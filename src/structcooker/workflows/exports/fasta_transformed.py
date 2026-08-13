from datacooker import RecipeBook

from structcooker.instructions.transforms.cifmol import convert_to_cifmol_transformed
from structcooker.instructions.transforms.sequence import build_fasta
from structcooker.mols import CIFMol

"""Build a fasta from a *transformed* CIF lmdb.

Like ``exports/fasta.py`` but the source records are already fanned into
``{cif_key: {cifmol_dict: ...}}`` (a rebuild output such as cif_valid_1), so they are
read with ``convert_to_cifmol_transformed`` rather than the base assembly-dict adapter.
"""

fasta_recipe = RecipeBook()

fasta_recipe.step(
    outputs=(("cifmol_dict", dict[str, dict[str, CIFMol]]),),
    instruction=convert_to_cifmol_transformed,
    kwargs={"value": ("db_data", dict)},
)

fasta_recipe.step(
    outputs=(("fasta", str),),
    instruction=build_fasta,
    kwargs={"cifmol_dict": ("cifmol_dict", dict[str, dict[str, CIFMol]] | None)},
)

RECIPE = fasta_recipe
TARGETS = ["fasta"]
