# CCD Ingest

Documentation for the **CCD ingest** workflow.

- Transform recipe: `structcooker.workflows.ingest.ccd` (`RECIPE`, `TARGETS`)
- Read / Write boundary: `configs/ingest/ccd_lmdb.yaml`

The recipe is the *Transform* layer only. The *Read* (loader) and *Write*
(serializer) boundaries come from the ingest config, not the recipe — both are
shown together in the pipeline diagram below.

## Pipeline (Read -> Transform -> Write)

Raw CCD `*.cif` files are loaded by `get_cif_data` into the three raw tables the
recipe consumes; the recipe produces `chem_comp_dict`; `to_bytes` serializes it
into the `biomol_CCD` LMDB.

![CCD ingest pipeline](sources/pipeline.png)

| Boundary | Hook | Value |
|----------|------|-------|
| Read  | `loader`     | `structcooker.instructions.readers.cif.get_cif_data` |
| Write | `serializer` | `structcooker.instructions.transforms.codecs.to_bytes` |
| Source | files   | `/public_data/CCD/components/*.cif*` |
| Dest   | LMDB    | `/public_data/CCD/biomol_CCD_202602.lmdb` |

## Transform recipe (detail)

![CCD recipe graph](sources/recipe.png)

```text
Targets: chem_comp_dict
Required inputs: <none>
Execution order:
1. _chem_comp_dict <- _group_rows(_chem_comp)
2. _chem_comp_atom_dict <- _group_rows(_chem_comp_atom)
3. _chem_comp_bond_dict <- _group_rows(_chem_comp_bond)
4. chem_comp_dict <- _parse_chem_comp(_chem_comp_dict, _chem_comp_atom_dict, _chem_comp_bond_dict)
```

### Mermaid (transform recipe)

```mermaid
flowchart LR
    classDef step fill:#f4f1e8,stroke:#6b5b3e,color:#2d2418;
    classDef target fill:#e6f0ff,stroke:#315ea8,color:#15233f;
    classDef input fill:#eef7ea,stroke:#4f7a39,color:#1d3112;
    classDef missing fill:#fdebec,stroke:#ba3d4c,color:#4a1218;
    classDef requested fill:#fff0c2,stroke:#b18100,color:#5e4300;
    input__chem_comp["_chem_comp"]
    input__chem_comp_atom["_chem_comp_atom"]
    input__chem_comp_bond["_chem_comp_bond"]
    step_1{{"step 1<br/>_chem_comp_dict <- _group_rows(_chem_comp)"}}
    step_2{{"step 2<br/>_chem_comp_atom_dict <- _group_rows(_chem_comp_atom)"}}
    step_3{{"step 3<br/>_chem_comp_bond_dict <- _group_rows(_chem_comp_bond)"}}
    step_4{{"step 4<br/>chem_comp_dict <- _parse_chem_comp(_chem_comp_dict, _chem_comp_atom_dict, _chem_comp_bond_dict)"}}
    target__chem_comp_dict["_chem_comp_dict"]
    target__chem_comp_atom_dict["_chem_comp_atom_dict"]
    target__chem_comp_bond_dict["_chem_comp_bond_dict"]
    target_chem_comp_dict["chem_comp_dict"]
    input__chem_comp --> step_1
    step_1 --> target__chem_comp_dict
    input__chem_comp_atom --> step_2
    step_2 --> target__chem_comp_atom_dict
    input__chem_comp_bond --> step_3
    step_3 --> target__chem_comp_bond_dict
    target__chem_comp_dict --> step_4
    target__chem_comp_atom_dict --> step_4
    target__chem_comp_bond_dict --> step_4
    step_4 --> target_chem_comp_dict
    class step_1 step;
    class step_2 step;
    class step_3 step;
    class step_4 step;
    class input__chem_comp input;
    class input__chem_comp_atom input;
    class input__chem_comp_bond input;
    class target__chem_comp_dict target;
    class target__chem_comp_atom_dict target;
    class target__chem_comp_bond_dict target;
    class target_chem_comp_dict requested;
```

## Sources

- `sources/pipeline.dot` / `.svg` / `.png` — Read -> Transform -> Write
- `sources/recipe.dot` / `.svg` / `.png` — Transform recipe DAG
