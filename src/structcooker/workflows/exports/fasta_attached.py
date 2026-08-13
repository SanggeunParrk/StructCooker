from datacooker import RecipeBook

from structcooker.instructions.transforms.cifmol import (
    convert_to_cifmol_attached_transformed,
)
from structcooker.instructions.transforms.sequence import build_fasta
from structcooker.mols import CIFMolAttached

"""Build a CIFMolAttached->fasta Cooker.

Same projection as ``exports/fasta.py`` but the source DB is an *attached* CIF lmdb
(schema B: each entry is ``{cif_key: {cifmol_attached_dict: ...}}``), so it reads the
records with ``convert_to_cifmol_attached_transformed`` instead of the base adapter.
Used to extract train/valid fastas that are header-consistent with this build's
``cif_pdb.fasta`` (production fastas carry different auth-chain labels).
"""

fasta_recipe = RecipeBook()

fasta_recipe.step(
    outputs=(("cifmol_dict", dict[str, dict[str, CIFMolAttached]]),),
    instruction=convert_to_cifmol_attached_transformed,
    kwargs={
        "value": ("db_data", dict),
    },
)

fasta_recipe.step(
    outputs=(("fasta", str),),
    instruction=build_fasta,
    kwargs={
        "cifmol_dict": ("cifmol_dict", dict[str, dict[str, CIFMolAttached]] | None),
    },
)

RECIPE = fasta_recipe
TARGETS = ["fasta"]
