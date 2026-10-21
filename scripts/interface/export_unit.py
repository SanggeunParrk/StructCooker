"""Write interface units as mmCIF for visual checks.

usage: export_unit.py <db_root> <out_dir> <unit_id> [<unit_id> ...] [--keep-bfactor]
In the output B_iso is a marker (unless --keep-bfactor): 100 interface residue, 50 other atom
within the contact cutoff of another node, 20 water (80 bridging), 0 the rest. In PyMOL:
``spectrum b, white_red`` or ``select iface, b > 90``. Waters are HOH with chain id W.
"""
import sys
from pathlib import Path

from structcooker.interface.db import InterfaceDB

ELEM = {"protein": "ATOM", "na": "ATOM"}


def to_cif(rec: dict, keep_b: bool = False) -> str:
    meta, nodes, at, rs, w = rec["meta"], rec["nodes"], rec["atoms"], rec["residues"], rec["water"]
    lines = [f"data_{meta['unit_id']}", f"# interface unit {meta['unit_id']}  res={meta['resolution']}  "
             f"nodes={'|'.join(nodes['kind'])}", "loop_", "_atom_site.group_PDB", "_atom_site.id",
             "_atom_site.type_symbol", "_atom_site.label_atom_id", "_atom_site.label_comp_id",
             "_atom_site.label_asym_id", "_atom_site.label_seq_id", "_atom_site.auth_seq_id",
             "_atom_site.auth_asym_id", "_atom_site.Cartn_x", "_atom_site.Cartn_y", "_atom_site.Cartn_z",
             "_atom_site.occupancy", "_atom_site.B_iso_or_equiv"]
    n = 0
    for i in range(len(at["xyz"])):
        r = int(at["res"][i])
        node = int(at["node"][i])
        kind = nodes["kind"][node]
        if keep_b:
            b = float(at["b_factor"][i])
        else:
            b = 100.0 if rs["interface"][r] else (50.0 if at["contact"][i] else 0.0)
        n += 1
        x, y, z = at["xyz"][i]
        lines.append(f"{ELEM.get(kind, 'HETATM')} {n} {at['element'][i]} {at['atom_name'][i].replace(chr(39), '_')} "
                     f"{rs['chem_comp_id'][r]} {nodes['chain_id'][node]} {rs['cif_idx'][r]} "
                     f"{rs['auth_idx'][r].split('.')[0]} {nodes['auth_asym_id'][node]} {x:.3f} {y:.3f} {z:.3f} "
                     f"{float(at['occupancy'][i]):.2f} {b:.2f}")
    for j in range(len(w["xyz"])):
        n += 1
        x, y, z = w["xyz"][j]
        b = float(80 if w["bridging"][j] else 20)
        lines.append(f"HETATM {n} O O HOH W {j + 1} {j + 1} W {x:.3f} {y:.3f} {z:.3f} 1.00 {b:.2f}")
    return "\n".join(lines) + "\n"


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    db, out, ids = InterfaceDB(args[0]), Path(args[1]), args[2:]
    out.mkdir(parents=True, exist_ok=True)
    for uid in ids:
        rec = db.get(uid)
        (out / f"{uid}.cif").write_text(to_cif(rec, "--keep-bfactor" in sys.argv))
        print(uid, rec["meta"]["unit_id"], len(rec["atoms"]["xyz"]), "atoms", len(rec["water"]["xyz"]), "waters")


if __name__ == "__main__":
    main()
