from pathlib import Path
from typing import TYPE_CHECKING, cast

import click
import lmdb
import networkx as nx
from datacooker.core import rebuild
from joblib import Parallel, delayed

from pipelines.cifmol import CIFMol, CIFMolAttached
from pipelines.utils import load_config
from pipelines.utils.convert import from_bytes
from pipelines.utils.lmdb import extract_key_list

if TYPE_CHECKING:
    from biomol.core.types import BioMolDict

_ENV_CACHE: dict[str, lmdb.Environment] = {}


def load_cif(key: str, env_path: Path) -> dict[str, dict[str, CIFMol]]:
    """
    Read a value from the LMDB database by key.

    Args:
        env_path: Path to the LMDB environment.
        key: Key of the data to retrieve.

    Returns
    -------
        dict
            The data dictionary retrieved from the LMDB database.
    """
    env_key = str(env_path)
    env = _ENV_CACHE.get(env_key)

    if env is None:
        env = lmdb.open(
            env_key,
            readonly=True,
            lock=False,
            max_readers=4096,
            readahead=True,
        )
        _ENV_CACHE[env_key] = env

    with env.begin(buffers=True) as txn:
        value = txn.get(key.encode())

    if value is None:
        msg = f"Key '{key}' not found in LMDB database at '{env_path}'."
        raise KeyError(msg)

    value = from_bytes(bytes(value))
    value, metadata = value["assembly_dict"], value["metadata_dict"]

    cifmol_dict: dict[str, dict[str, CIFMol]] = {}
    for cif_key, _item in value.items():
        assembly_id, model_id, alt_id = cif_key.split("_")

        md = dict(metadata)
        md["assembly_id"] = assembly_id
        md["model_id"] = model_id
        md["alt_id"] = alt_id

        item = dict(_item)
        item["metadata"] = md
        item = cast("BioMolDict", item)

        cifmol_dict[cif_key] = {"cifmol": CIFMol.from_dict(item)}

    return cifmol_dict


def load_cifmol_attached(
    key: str,
    env_path: Path,
) -> dict[str, dict[str, CIFMolAttached]]:
    """Read a CIFMolAttached object from the LMDB database by cifID."""
    env_key = str(env_path)
    env = _ENV_CACHE.get(env_key)
    if env is None:
        env = lmdb.open(
            env_key,
            readonly=True,
            lock=False,
            max_readers=4096,
            readahead=True,
        )
        _ENV_CACHE[env_key] = env

    with env.begin(buffers=True) as txn:
        value = txn.get(key.encode())

    if value is None:
        msg = f"Key '{key}' not found in LMDB database at '{env_path}'."
        raise KeyError(msg)

    value = from_bytes(bytes(value))

    cifmol_dict: dict[str, dict[str, CIFMolAttached]] = {}
    for cif_key, _item in value.items():
        if _item["cifmol_dict"] is None:
            continue
        cifmol_dict[cif_key] = {
            "cifmol": CIFMolAttached.from_dict(_item["cifmol_dict"]),
        }
    return cifmol_dict


def load_fasta(fasta_path: Path) -> dict[str, str]:
    """Load fasta file into a dictionary."""
    fasta_dict = {}
    with fasta_path.open("r") as f:
        lines = f.readlines()
    current_header = ""
    for _line in lines:
        line = _line.strip()
        if line.startswith(">"):
            current_header = line[1:]
            fasta_dict[current_header] = ""
        else:
            fasta_dict[current_header] += line
    return fasta_dict


def load_a3m(a3m_path: Path) -> dict[str, str]:
    """Load a3m file into a dictionary."""
    a3m_dict = {}
    with a3m_path.open("r") as f:
        lines = f.readlines()
    current_header = ""
    for _line in lines:
        line = _line.strip()
        if line.startswith(">"):
            current_header = line[1:]
            a3m_dict[current_header] = ""
        else:
            a3m_dict[current_header] += line
    return a3m_dict


def load_graph(key: str, env_path: Path) -> nx.Graph:
    """
    Read a value from the LMDB database by key.

    Args:
        env_path: Path to the LMDB environment.
        key: Key of the data to retrieve.

    Returns
    -------
        dict
            The data dictionary retrieved from the LMDB database.
    """
    env_key = str(env_path)

    env = _ENV_CACHE.get(env_key)
    if env is None:
        env = lmdb.open(
            env_key,
            readonly=True,
            lock=False,
            max_readers=4096,
            readahead=True,
        )
        _ENV_CACHE[env_key] = env

    with env.begin(buffers=True) as txn:
        value = txn.get(key.encode())

    if value is None:
        msg = f"Key '{key}' not found in LMDB database at '{env_path}'."
        raise KeyError(msg)

    value = from_bytes(bytes(value))["cluster_graph"]
    chain_id_list = value["nodes"]["chain_ids"]["value"]
    seq_clusters = value["nodes"]["seq_clusters"]["value"]
    edge_list = value["edges"]["contact_edges"]
    src, dst = edge_list["src_indices"], edge_list["dst_indices"]
    graph = nx.Graph()

    for chain_id, cluster_id in zip(chain_id_list, seq_clusters, strict=True):
        graph.add_node(chain_id, label=cluster_id)

    for s, d in zip(src, dst, strict=True):
        graph.add_edge(chain_id_list[s], chain_id_list[d])

    return graph


def cifmol_process(
    cif_id: str,
    recipe_path: Path,
    cif_db_path: Path,
    targets: list[str] | None = None,
) -> dict:
    """Parse a CIF file using a predefined recipe."""
    cifmol_dict = load_cif(cif_id, env_path=cif_db_path)
    output = {}
    for key in cifmol_dict:
        inner_cifmol_dict = cifmol_dict[key]
        output[key] = rebuild(
            datadict=inner_cifmol_dict,
            recipe_path=recipe_path,
            targets=targets,
        )
    return output


@click.group()
def cli() -> None:
    """Build and merge LMDB databases from CIF files."""


@cli.command("extract_fasta")
@click.argument("config", type=click.Path(exists=True, path_type=Path))
@click.option("--njobs", "-j", type=int, default=-1, show_default=True)
def extract_fasta(
    config: Path,
    njobs: int = -1,
) -> None:
    """Extract fasta from CIFMol objects stored in LMDB database."""
    # python -m scripts.postprocess extract_fasta ~/data/BioMolDBv2_2024Oct21/cif.lmdb/ ~/data/BioMolDBv2_2024Oct21/fasta/
    config_dict = load_config(config)
    key_list = extract_key_list(cif_db_path)
    click.echo(f"Extracting fasta for {len(key_list)} CIFMol objects...")
    from pipelines.recipe.extract_fasta import RECIPE, TARGETS

    recipe, targets = RECIPE, TARGETS

    results = Parallel(n_jobs=njobs, verbose=10)(
        delayed(cifmol_process)(
            cif_id=key,
            recipe=recipe,
            targets=targets,
            cif_db_path=cif_db_path,
        )
        for key in key_list
    )
    results = dict(zip(key_list, results, strict=False))

    # remove redundancy
    path_to_lines = {}
    all_fasta = []

    for cif_key, outer_data in results.items():
        if outer_data is None:
            continue
        for inner_data in outer_data.values():
            for chain_id in inner_data["fasta"]:
                path = output_path / cif_key[1:3] / f"{cif_key}_{chain_id}.fasta"
                path_to_lines[path] = inner_data["fasta"][chain_id]
                all_fasta.append(inner_data["fasta"][chain_id])

    def _write_fasta(
        output_path: Path,
        lines: str,
    ) -> None:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with output_path.open("w") as f:
            f.writelines(lines)

    click.echo(f"Writing {len(path_to_lines)} fasta files to {output_path}...")
    Parallel(n_jobs=njobs, verbose=10)(
        delayed(_write_fasta)(
            output_path=path,
            lines=lines,
        )
        for path, lines in path_to_lines.items()
    )

    # write merged_fasta file
    merged_fasta_path = output_path / "merged.fasta"
    click.echo(f"Writing merged fasta file to {merged_fasta_path}...")
    _write_fasta(
        output_path=merged_fasta_path,
        lines="".join(all_fasta),
    )


@cli.command("build_seq_hash_map")
@click.argument("merged_fasta_path", type=click.Path(path_type=Path))
@click.argument("output_path", type=click.Path(path_type=Path))
def build_seq_hash_map(
    merged_fasta_path: Path,
    output_path: Path,
) -> None:
    """Build a sequence hash map from the merged fasta file."""
    # python -m scripts.postprocessing build_seq_hash_map /public_data/BioMolDBv2_2024Oct21/fasta/merged.fasta /public_data/BioMolDBv2_2024Oct21/metadata/seq_hash_map.tsv
    fasta_dict = load_fasta(merged_fasta_path)
    results = rebuild(
        datadict={"fasta_dict": fasta_dict},
        recipe_path=recipe_path,
    )
    seq_hash_map = results["seq_hash_map"]

    # sort by 1. identifier, 2. sequence length
    seq_hash_map = dict(
        sorted(
            seq_hash_map.items(),
            key=lambda item: (item[0][0], len(item[1])),
        ),
    )

    # write to output_path
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w") as f:
        for seq_hash, sequence in seq_hash_map.items():
            f.write(f"{seq_hash}\t{sequence}\n")


@cli.command("seq_cluster")
@click.argument("seq_hash_map", type=click.Path(path_type=Path))
@click.argument("sabdab_summary_path", type=click.Path(path_type=Path))
@click.argument("tmp_dir", type=click.Path(path_type=Path))
@click.argument("output_path", type=click.Path(path_type=Path))
@click.option("--mmseqs2_seq_id", type=float, default=0.3, show_default=True)
@click.option("--mmseqs2_cov", type=float, default=0.8, show_default=True)
@click.option("--mmseqs2_covmode", type=int, default=0, show_default=True)
@click.option("--mmseqs2_clustermode", type=int, default=1, show_default=True)
def seq_cluster(  # noqa: PLR0913
    seq_hash_map: Path,
    sabdab_summary_path: Path,
    tmp_dir: Path,
    output_path: Path,
    mmseqs2_seq_id: float,
    mmseqs2_cov: float,
    mmseqs2_covmode: int,
    mmseqs2_clustermode: int,
) -> None:
    """
    Cluster sequences...

    - Antibody: SabDab, CD-HIT (H3 if exists, else L3) cluster ID = A1234567
    - Peptides <10 residues, 100%
    - Polypeptide(D) (따로 나눠서)
    - Other Proteins: MMseqs2 easy-cluster (seq_id=0.3,cov=0.8,covmode=0,clustmode=1) cluster ID = P1234567
    - RNA, DNA, Ligands: 100% sequence identity cluster ID = N1234567 (same as sequence hash)
        /public_data/BioMolDBv2_2024Oct21/metadata/seq_hash_map.tsv ./postprocessing/tmp/ /public_data/BioMolDBv2_2024Oct21/cluster/seq_clusters.tsv

    Example:
    python -m scripts.postprocessing seq_cluster \
        /public_data/BioMolDBv2_2024Oct21/metadata/seq_hash_map.tsv \
        ./postprocessing/external_source/SabDab/sabdab_summary_all.tsv \
        /public_data/BioMolDBv2_2024Oct21/seq_cluster/tmp/ \
        /public_data/BioMolDBv2_2024Oct21/seq_cluster/seq_clusters.tsv \
        --mmseqs2_seq_id 0.3 \
        --mmseqs2_cov 0.8 \
        --mmseqs2_covmode 0 \
        --mmseqs2_clustermode 1
    """
    # tmp dir
    tmp_dir.mkdir(parents=True, exist_ok=True)
    # recipe path
    from pipelines.recipe.seq_cluster import RECIPE, TARGETS

    recipe, targets = RECIPE, TARGETS
    cluster_dict = rebuild(
        datadict={
            "tmp_dir": tmp_dir,
            "seq_hash_map": seq_hash_map,
            "sabdab_summary_path": sabdab_summary_path,
            "mmseqs2_params": {
                "mmseqs2_seq_id": mmseqs2_seq_id,
                "mmseqs2_cov": mmseqs2_cov,
                "mmseqs2_covmode": str(mmseqs2_covmode),
                "mmseqs2_clustermode": str(mmseqs2_clustermode),
            },
        },
        recipe_path=recipe_path,
    )

    # write down the cluster_dict
    output_path.parent.mkdir(parents=True, exist_ok=True)
    cluster_num = len(cluster_dict["cluster_dict"])
    total_members = 0
    with output_path.open("w") as f:
        for rep_seq_hash, member_list in cluster_dict["cluster_dict"].items():
            members = ",".join(member_list)
            f.write(f"c{rep_seq_hash}\t{members}\n")
            total_members += len(member_list)
    click.echo(f"Saved {cluster_num} clusters with {total_members} total members.")


@cli.command("graph_lmdb_attached")
@click.argument("cif_db_path", type=click.Path(path_type=Path))
@click.argument("graph_lmdb_path", type=click.Path(path_type=Path))
def graph_lmdb_attached(
    cif_db_path: Path,
    graph_lmdb_path: Path,
) -> None:
    """
    Extract graphs of CIFMol objects stored in LMDB database.

    Example:
    python -m scripts.postprocessing graph_lmdb \
        /public_data/BioMolDBv2_2024Oct21/cif.lmdb \
        /public_data/BioMolDBv2_2024Oct21/metadata/seq_hash_map.tsv \
        /public_data/BioMolDBv2_2024Oct21/seq_cluster/seq_clusters.tsv \
        /public_data/BioMolDBv2_2024Oct21/graph.lmdb
    """
    key_list = extract_key_list(cif_db_path)
    from pipelines.recipe.graph_lmdb_from_attached import RECIPE, TARGETS

    recipe, targets = RECIPE, TARGETS

    # ---------------------------------------
    # Worker on CHUNK basis
    # ---------------------------------------
    def _process_chunk(keys: list[str]) -> dict[str, bytes]:
        out: dict[str, bytes] = {}
        for cif_id in keys:
            cifmol_dict = load_cifmol_attached(cif_id, env_path=cif_db_path)
            if len(cifmol_dict) == 0:
                continue
            for inner_key, obj in cifmol_dict.items():
                cifmol = obj["cifmol"]
                result = rebuild(
                    datadict={
                        "cifmol": cifmol,
                    },
                    recipe_path=recipe_path,
                )
                graph_bytes = result["graph_bytes"]
                full_key = f"{cif_id}_{inner_key}"
                out[full_key] = graph_bytes
        return out

    # make chunks
    CHUNK = 200  # tuneable
    key_chunks = [key_list[i : i + CHUNK] for i in range(0, len(key_list), CHUNK)]
    click.echo(f"Processing {len(key_chunks)} chunks of size {CHUNK}...")

    # parallel batch
    chunk_results = Parallel(n_jobs=-1, verbose=10)(
        delayed(_process_chunk)(chunk) for chunk in key_chunks
    )

    # merge
    graph_dict: dict[str, bytes] = {}
    for d in chunk_results:
        if d is None:
            continue
        graph_dict.update(d)

    # write to lmdb (original behaviour)
    all_keys = list(graph_dict.keys())
    click.echo(f"Writing {len(graph_dict)} graphs to LMDB at {graph_lmdb_path}...")
    env = lmdb.open(str(graph_lmdb_path), map_size=int(1e12))

    WRITE_CHUNK = 10_000
    for i in range(0, len(all_keys), WRITE_CHUNK):
        click.echo(
            f"Processing files {i} to {min(i + WRITE_CHUNK, len(all_keys))} / {len(all_keys)}",
        )
        data_chunk = all_keys[i : i + WRITE_CHUNK]
        with env.begin(write=True) as txn:
            for key in data_chunk:
                zcompressed = graph_dict[key]
                txn.put(key.encode(), zcompressed)


@cli.command("extract_edge")
@click.argument("graph_db_path", type=click.Path(path_type=Path))
@click.argument("output_path", type=click.Path(path_type=Path))  # tsv file
def extract_edge(
    graph_db_path: Path,
    output_path: Path,
) -> None:
    """
    Extract edges from graph LMDB database.

    Example:
    python -m scripts.postprocess extract_edge \
        ~/data/BioMolDBv2_2024Oct21/graph.lmdb \
        ~/data/BioMolDBv2_2024Oct21/graph_edges.tsv
    """
    key_list = extract_key_list(graph_db_path)

    # ---------------------------------------
    # Worker on CHUNK basis
    # ---------------------------------------
    def _process_chunk(keys: list[str]) -> dict[tuple[str, str], list[str]]:
        edges_dict: dict[tuple[str, str], list[str]] = {}
        for cif_id in keys:
            graph = load_graph(cif_id, env_path=graph_db_path)

            for c1, c2 in graph.edges():
                sc1 = graph.nodes[c1]["label"]
                sc2 = graph.nodes[c2]["label"]

                if sc1 <= sc2:
                    key = (sc1, sc2)
                    id_str = f"{cif_id}_({c1})_({c2})"
                else:
                    key = (sc2, sc1)
                    id_str = f"{cif_id}_({c2})_({c1})"
                if key not in edges_dict:
                    edges_dict[key] = []
                edges_dict[key].append(id_str)

        return edges_dict

    # make chunks
    CHUNK = 200  # tuneable
    key_chunks = [key_list[i : i + CHUNK] for i in range(0, len(key_list), CHUNK)]
    click.echo(f"Processing {len(key_chunks)} chunks of size {CHUNK}...")

    # parallel batch
    chunk_results = Parallel(n_jobs=-1, verbose=10)(
        delayed(_process_chunk)(chunk) for chunk in key_chunks
    )

    # merge
    merged_edges_dict: dict[tuple[str, str], list[str]] = {}
    for d in chunk_results:
        if d is None:
            continue
        for key, id_list in d.items():
            if key not in merged_edges_dict:
                merged_edges_dict[key] = []
            merged_edges_dict[key].extend(id_list)

    # write to output_path
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w") as f:
        for (c1, c2), id_list in merged_edges_dict.items():
            ids = ",".join(id_list)
            f.write(f"{c1}\t{c2}\t{ids}\n")


if __name__ == "__main__":
    cli()
