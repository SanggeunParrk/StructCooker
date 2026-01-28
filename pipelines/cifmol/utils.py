from pathlib import Path

import numpy as np

from pipelines.cifmol import CIFMol, CIFMolAttached


def to_cif(cifmol: CIFMol|CIFMolAttached, output_path: Path) -> None:  # noqa: PLR0915
    """
    Write a CIFMol object to a CIF file.

    1. _atom_site
    2. _struct_conn (TODO)
    3. ~_scheme (TODO)
    4. ~_branch (TODO)
    5. ~_entity (TODO)
    """

    def _to_mmcif_format(array: np.ndarray) -> list:
        _list = [str(item) for item in array]
        max_length = max([len(item) for item in _list])
        return [item.ljust(max_length) for item in _list]

    output = f"#\ndata_{cifmol.id}_{cifmol.model_id}_{cifmol.alt_id}\n"

    header = [
        "_atom_site.group_PDB",
        "_atom_site.id",
        "_atom_site.type_symbol",
        "_atom_site.label_atom_id",
        "_atom_site.label_alt_id",
        "_atom_site.label_comp_id",
        "_atom_site.label_asym_id",
        "_atom_site.label_entity_id",
        "_atom_site.label_seq_id",
        "_atom_site.pdbx_PDB_ins_code",
        "_atom_site.Cartn_x",
        "_atom_site.Cartn_y",
        "_atom_site.Cartn_z",
        "_atom_site.occupancy",
        "_atom_site.B_iso_or_equiv",
        "_atom_site.pdbx_formal_charge",
        "_atom_site.auth_seq_id",
        "_atom_site.auth_comp_id",
        "_atom_site.auth_asym_id",
        "_atom_site.auth_atom_id",
        "_atom_site.pdbx_PDB_model_num",
    ]
    output += "#\n"
    output += "loop_\n"
    output += "\n".join(header) + "\n"

    xyz = cifmol.atoms.xyz
    mask = ~np.isnan(xyz).value.any(axis=1)
    length = mask.sum()
    atom_to_res = np.array(cifmol.index_table.atom_to_res)
    res_to_chain = np.array(cifmol.index_table.res_to_chain)
    atom_to_chain = res_to_chain[atom_to_res]

    group_PDB_list = cifmol.residues.hetero[atom_to_res][mask].value  # noqa: N806
    id_list = 1 + np.arange(length)
    type_symbol_list = cifmol.atoms.element[mask].value
    # label_atom_id_list = cifmol.atoms.id[mask].value
    label_atom_id_list = cifmol.atoms.id[mask].value
    label_alt_id_list = [cifmol.alt_id] * length
    label_comp_id_list = cifmol.residues.chem_comp_id[atom_to_res][mask].value
    label_asym_id_list = cifmol.chains.chain_id[atom_to_chain][mask].value
    label_entity_id_list = cifmol.chains.entity_id[atom_to_chain][mask].value
    label_seq_id_list = cifmol.residues.cif_idx[atom_to_res][mask].value
    auth_idx_list = cifmol.residues.auth_idx[atom_to_res][mask].value
    auth_seq_id_list, ins_code_list = zip(
        *[
            (p[0], p[1] if len(p) > 1 else "?")
            for p in (s.split(".") for s in auth_idx_list)
        ],
        strict=False,
    )

    cartn_x_list = cifmol.atoms.xyz[mask, 0].value
    cartn_y_list = cifmol.atoms.xyz[mask, 1].value
    cartn_z_list = cifmol.atoms.xyz[mask, 2].value
    occupancy_list = cifmol.atoms.occupancy[mask].value
    b_iso_or_equiv_list = cifmol.atoms.b_factor[mask].value
    pdbx_formal_charge_list = cifmol.atoms.charge[mask].value
    pdbx_PDB_model_num_list = [cifmol.model_id] * length  # noqa: N806

    # to mmcif format
    group_PDB_list = _to_mmcif_format(group_PDB_list)  # noqa: N806
    type_symbol_list = _to_mmcif_format(type_symbol_list)
    label_atom_id_list = _to_mmcif_format(label_atom_id_list)
    label_comp_id_list = _to_mmcif_format(label_comp_id_list)
    label_asym_id_list = _to_mmcif_format(label_asym_id_list)
    label_entity_id_list = _to_mmcif_format(label_entity_id_list)
    label_seq_id_list = _to_mmcif_format(label_seq_id_list)
    ins_code_list = _to_mmcif_format(ins_code_list)
    auth_idx_list = _to_mmcif_format(auth_idx_list)
    auth_seq_id_list = _to_mmcif_format(auth_seq_id_list)
    id_list = _to_mmcif_format(id_list)
    cartn_x_list = [f"{x:.3f}" for x in cartn_x_list]
    cartn_y_list = [f"{y:.3f}" for y in cartn_y_list]
    cartn_z_list = [f"{z:.3f}" for z in cartn_z_list]
    cartn_x_list = _to_mmcif_format(cartn_x_list)
    cartn_y_list = _to_mmcif_format(cartn_y_list)
    cartn_z_list = _to_mmcif_format(cartn_z_list)
    occupancy_list = [f"{o:.2f}" for o in occupancy_list]
    b_iso_or_equiv_list = [f"{b:.2f}" for b in b_iso_or_equiv_list]
    occupancy_list = _to_mmcif_format(occupancy_list)
    b_iso_or_equiv_list = _to_mmcif_format(b_iso_or_equiv_list)
    pdbx_formal_charge_list = _to_mmcif_format(pdbx_formal_charge_list)
    pdbx_PDB_model_num_list = _to_mmcif_format(pdbx_PDB_model_num_list)  # noqa: N806

    for idx in range(length):
        fields = [
            group_PDB_list[idx],
            id_list[idx],
            type_symbol_list[idx],
            label_atom_id_list[idx],
            label_alt_id_list[idx],
            label_comp_id_list[idx],
            label_asym_id_list[idx],
            label_entity_id_list[idx],
            label_seq_id_list[idx],
            ins_code_list[idx],
            cartn_x_list[idx],
            cartn_y_list[idx],
            cartn_z_list[idx],
            occupancy_list[idx],
            b_iso_or_equiv_list[idx],
            pdbx_formal_charge_list[idx],
            auth_seq_id_list[idx],
            label_comp_id_list[idx],
            label_asym_id_list[idx],
            label_atom_id_list[idx],
            pdbx_PDB_model_num_list[idx],
        ]
        output += " ".join(map(str, fields)) + "\n"

    output_path.write_text(output)
