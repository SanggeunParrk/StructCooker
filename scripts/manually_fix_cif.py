from pathlib import Path
import click
import os

@click.group()
def cli():
	pass

@cli.command("manually-fix")
@click.argument('src_dir', type=click.Path(exists=True, file_okay=False))
@click.argument('error_item_dir', type=click.Path(exists=True, file_okay=False))
@click.argument('dst_dir', type=click.Path(exists=True, file_okay=False))
def manually_fix(src_dir: Path, error_item_dir: Path, dst_dir: Path):
	"""
	Manually copy files from src to dst, skipping files that already exist.
	"""
	src_dir = Path(src_dir)
	dst_dir = Path(dst_dir)
	error_item_dir = Path(error_item_dir)
	file_list = list(error_item_dir.glob('*.cif.gz'))
	pdb_list = [file.stem[:4] for file in file_list]
	for pdb_id in pdb_list:
		file = src_dir / f"{pdb_id[1:3]}" / f"{pdb_id}.cif.gz"
		inner_dir = dst_dir/ Path(pdb_id[1:3])
		dst_file = inner_dir / f"{pdb_id}.cif.gz"
		breakpoint()
		if not dst_file.exists():
			breakpoint()
		if not inner_dir.exists():
			inner_dir.mkdir(parents=True, exist_ok=True)
		os.system(f'cp {file} {dst_file}')

if __name__ == '__main__':
	cli()
	# python manually_fix_cif.py manually-fix /public_data/BioMolDB_2024Oct21/cif/cif_raw/ /home/psk6950/data/BioMolDB_20260224/cif/error_items/manually_fixed /home/psk6950/data/BioMolDB_20260224/cif/raw/