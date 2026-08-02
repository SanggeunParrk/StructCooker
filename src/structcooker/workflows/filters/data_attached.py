from datetime import date

from datacooker import RecipeBook

from structcooker.instructions.transforms.filtering import (
    filter_by_resolution_and_date,
    filter_cifmol_by_polymer_chain_count,
    filter_cifmol_by_token_count,
    filter_signalp,
)
from structcooker.instructions.transforms.sequence import filter_water
from structcooker.mols import CIFMol

"""Train-filter an ALREADY-ATTACHED CIF lmdb (AF3 training subset).

Same filtering pipeline as ``filters/data.py`` (token/chain/resolution-date/
water/signalp trim) but the input is a ``CIFMolAttached`` (seq_id + cluster_id
already assigned on the full, untrimmed sequences) and the final output target
is ``cifmol_attached_dict`` so the trimmed records keep their cluster_id and can
be read back with ``convert_to_cifmol_attached_transformed``.

Attaching on full sequences and filtering afterwards is the correct order:
signalp trimming shortens chain sequences, so attaching *after* the trim breaks
the seq->seq_id lookup (the trimmed sequence is not a key in seq_id_map). The
attach in cif_pdb_attached.lmdb was done on the untrimmed sequences; the signalp
lookup here still reads the untrimmed chain sequences, so it resolves fine.
"""

recipe = RecipeBook()


recipe.step(
    outputs=(("cifmol_filtered_by_token_count", CIFMol),),
    instruction=filter_cifmol_by_token_count,
    kwargs={
        "cifmol": ("cifmol", CIFMol),
        "min_token_count": ("min_token_count", int),
        "max_token_count": ("max_token_count", int | None),
    },
)


recipe.step(
    outputs=(("cifmol_filtered_by_chain_count", CIFMol),),
    instruction=filter_cifmol_by_polymer_chain_count,
    kwargs={
        "cifmol": ("cifmol_filtered_by_token_count", CIFMol),
        "max_polymer_chain_count": ("max_polymer_chain_count", int),
    },
)


recipe.step(
    outputs=(("cifmol_filtered_by_resolution_date", CIFMol),),
    instruction=filter_by_resolution_and_date,
    kwargs={
        "resolution_cutoff": ("resolution_cutoff", float),
        "start_date": ("start_date", date | str),
        "end_date": ("end_date", date | str),
        "cifmol": ("cifmol_filtered_by_chain_count", CIFMol),
    },
)


recipe.step(
    outputs=(("cifmol_wo_water", CIFMol),),
    instruction=filter_water,
    kwargs={
        "cifmol": ("cifmol_filtered_by_resolution_date", CIFMol),
    },
)


recipe.step(
    outputs=(("cifmol_attached_dict", dict),),
    instruction=filter_signalp,
    kwargs={
        "cifmol": ("cifmol_wo_water", CIFMol),
        "seqid_map": ("seqid_map", dict),
        "signalp_dict": ("signalp_dict", dict),
    },
)


RECIPE = recipe
TARGETS = ["cifmol_attached_dict"]
