# ruff: noqa: E402

import ast
import importlib.util
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import click
import lmdb
from datacooker import parse_dict
from omegaconf import OmegaConf

from pipelines.utils import load_config
from pipelines.utils.lmdb import read_lmdb

OmegaConf.register_new_resolver("p", lambda x: Path(x))


@dataclass(frozen=True)
class RecipeTestConfig:
    """Runtime options for recipe validation."""

    db_path: Path
    recipe_path: Path
    db_field: str
    targets: list[str] | None
    config_dict: dict[str, Any]
    extra_inputs: dict[str, Any]


def _load_recipe_targets(recipe_path: Path) -> list[str] | None:
    """Load TARGETS from a recipe module if present."""
    module_name = f"test_db_{recipe_path.stem}"
    spec = importlib.util.spec_from_file_location(module_name, recipe_path)
    if spec is None or spec.loader is None:
        msg = f"Failed to load recipe module from '{recipe_path}'."
        raise ImportError(msg)

    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    targets = getattr(module, "TARGETS", None)
    if targets is None:
        return None
    return list(targets)


def _parse_cli_value(raw_value: str) -> object:
    """Parse CLI input values into Python objects when possible."""
    try:
        return ast.literal_eval(raw_value)
    except (SyntaxError, ValueError):
        if raw_value.startswith(("/", "./", "../")) or "/" in raw_value:
            return Path(raw_value)
        return raw_value


def _parse_assignments(items: tuple[str, ...]) -> dict[str, Any]:
    """Parse repeated KEY=VALUE CLI arguments."""
    parsed: dict[str, Any] = {}
    for item in items:
        if "=" not in item:
            msg = f"Invalid assignment '{item}'. Expected KEY=VALUE."
            raise click.ClickException(msg)
        key, raw_value = item.split("=", 1)
        if not key:
            msg = f"Invalid assignment '{item}'. KEY must not be empty."
            raise click.ClickException(msg)
        parsed[key] = _parse_cli_value(raw_value)
    return parsed


def _get_entry_count(db_path: Path) -> int:
    """Return the number of entries in the LMDB."""
    env = lmdb.open(str(db_path), readonly=True, lock=False)
    try:
        return env.stat()["entries"]
    finally:
        env.close()


def _get_sample_keys(
    db_path: Path,
    sample_size: int,
) -> list[str]:
    """Collect the first N keys from the LMDB for a lightweight test."""
    env = lmdb.open(str(db_path), readonly=True, lock=False)
    try:
        keys: list[str] = []
        with env.begin() as txn:
            cursor = txn.cursor()
            for index, key in enumerate(
                cursor.iternext(keys=True, values=False),
                start=1,
            ):
                keys.append(key.decode())
                if index >= sample_size:
                    break
        return keys
    finally:
        env.close()


def _resolve_db_path(config_dict: dict[str, Any], cli_db_path: Path | None) -> Path:
    """Resolve the LMDB path from CLI or config."""
    db_path = cli_db_path or config_dict.get("db_path") or config_dict.get("env_path")
    if db_path is None:
        msg = "DB path is required. Use --db-path or provide a config with db_path/env_path."
        raise click.ClickException(msg)

    db_path = Path(db_path)
    if not db_path.exists():
        msg = f"DB path does not exist: {db_path}"
        raise click.ClickException(msg)
    return db_path


def _resolve_recipe_path(
    config_dict: dict[str, Any],
    cli_recipe_path: Path | None,
) -> Path | None:
    """Resolve a recipe path intended for LMDB post-processing tests."""
    if cli_recipe_path is not None:
        return cli_recipe_path

    extract_recipe_path = config_dict.get("extract_recipe_path")
    if extract_recipe_path is not None:
        return Path(extract_recipe_path)

    recipe_path = config_dict.get("recipe_path")
    if recipe_path is not None and "db_path" in config_dict:
        return Path(recipe_path)

    return None


def _resolve_expected_keys(
    config_dict: dict[str, Any],
    cli_expected_keys: tuple[str, ...],
) -> list[str] | None:
    """Resolve expected top-level keys for structure validation."""
    if cli_expected_keys:
        return list(cli_expected_keys)

    build_recipe_path = config_dict.get("recipe")
    env_path = config_dict.get("env_path")
    if build_recipe_path is None or env_path is None:
        return None

    return _load_recipe_targets(Path(build_recipe_path))


def _print_header(title: str) -> None:
    click.echo(f"\n[{title}]")


def _run_structure_check(
    db_path: Path,
    keys: list[str],
    expected_top_keys: list[str] | None,
) -> int:
    """Check that sampled LMDB entries are readable and structurally consistent."""
    failures = 0
    signatures: Counter[tuple[str, ...]] = Counter()
    expected_set = set(expected_top_keys or [])

    _print_header("DB Structure")
    for key in keys:
        try:
            value = read_lmdb(db_path, key)
            top_keys = sorted(value.keys())
            signatures.update([tuple(top_keys)])

            missing_keys = sorted(expected_set - set(top_keys))
            if missing_keys:
                click.echo(
                    f"FAIL structure key={key} missing_keys={missing_keys} top_keys={top_keys}",
                )
                failures += 1
                continue

            click.echo(f"PASS structure key={key} top_keys={top_keys}")
        except Exception as error:  # noqa: BLE001
            click.echo(f"FAIL structure key={key} error={error}")
            failures += 1

    if len(signatures) > 1:
        click.echo(f"FAIL structure signatures={dict(signatures)}")
        failures += 1
    elif signatures:
        click.echo(f"PASS structure signature={next(iter(signatures))}")

    return failures


def _build_recipe_datadict(
    key: str,
    db_path: Path,
    db_field: str,
    config_dict: dict[str, Any],
    extra_inputs: dict[str, Any],
) -> dict[str, Any]:
    """Load one LMDB entry and package inputs for recipe execution."""
    data = read_lmdb(db_path, key)
    convert_func = config_dict.get("convert_func")
    transform_func = config_dict.get("transform_func")
    if convert_func is not None:
        data = convert_func(data)

    datadict = {db_field: data}
    datadict.update(config_dict)
    datadict.update(config_dict.get("additional_inputs", {}))
    datadict.update(extra_inputs)
    datadict["transform_func"] = transform_func
    return datadict


def _run_recipe_for_key(key: str, recipe_config: RecipeTestConfig) -> list[str]:
    """Execute one recipe run and return served result keys."""
    datadict = _build_recipe_datadict(
        key=key,
        db_path=recipe_config.db_path,
        db_field=recipe_config.db_field,
        config_dict=recipe_config.config_dict,
        extra_inputs=recipe_config.extra_inputs,
    )
    transform_func = recipe_config.config_dict.get("transform_func")
    result = parse_dict(
        recipe_path=recipe_config.recipe_path,
        datadict=datadict,
        transform_func=transform_func,
        targets=recipe_config.targets,
    )
    return sorted(result.keys())


def _collect_recipe_result(
    key: str,
    recipe_config: RecipeTestConfig,
) -> tuple[bool, list[str] | str]:
    """Execute one recipe run and collect either result keys or an error."""
    try:
        return True, _run_recipe_for_key(key, recipe_config)
    except Exception as error:  # noqa: BLE001
        return False, str(error)


def _run_recipe_check(keys: list[str], recipe_config: RecipeTestConfig) -> int:
    """Execute the recipe on sampled keys and print PASS/FAIL."""
    failures = 0
    resolved_targets = recipe_config.targets or _load_recipe_targets(
        recipe_config.recipe_path,
    )
    if resolved_targets is None:
        msg = (
            f"Recipe '{recipe_config.recipe_path}' does not declare TARGETS. "
            "Use --target to specify outputs."
        )
        raise click.ClickException(msg)
    recipe_config = RecipeTestConfig(
        db_path=recipe_config.db_path,
        recipe_path=recipe_config.recipe_path,
        db_field=recipe_config.db_field,
        targets=resolved_targets,
        config_dict=recipe_config.config_dict,
        extra_inputs=recipe_config.extra_inputs,
    )

    _print_header("Recipe Test")
    click.echo(f"recipe={recipe_config.recipe_path}")
    click.echo(f"targets={resolved_targets}")

    for key in keys:
        ok, payload = _collect_recipe_result(key, recipe_config)
        if ok:
            result_keys = payload
            click.echo(f"PASS recipe key={key} result_keys={result_keys}")
            continue

        click.echo(f"FAIL recipe key={key} error={payload}")
        failures += 1

    return failures


@click.command()
@click.option(
    "--config",
    "config_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="YAML config path. Supports existing configs/* structure in this repo.",
)
@click.option(
    "--db-path",
    type=click.Path(exists=True, path_type=Path),
    help="LMDB path to test. Overrides db_path/env_path from --config.",
)
@click.option(
    "--recipe-path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Recipe to run against sampled LMDB entries.",
)
@click.option(
    "--key",
    "keys",
    multiple=True,
    help="Specific LMDB key(s) to test. If omitted, the script samples the first N keys.",
)
@click.option(
    "--sample-size",
    type=click.IntRange(min=1),
    default=5,
    show_default=True,
    help="Number of keys to sample when --key is not provided.",
)
@click.option(
    "--expected-key",
    "expected_keys",
    multiple=True,
    help="Expected top-level keys in each LMDB entry.",
)
@click.option(
    "--target",
    "targets",
    multiple=True,
    help="Recipe target(s) to serve. If omitted, TARGETS from the recipe are used.",
)
@click.option(
    "--input",
    "input_items",
    multiple=True,
    help="Extra recipe inputs as KEY=VALUE.",
)
@click.option(
    "--db-field",
    default="db_data",
    show_default=True,
    help="Field name used when passing LMDB entries into the recipe.",
)
@click.option(
    "--skip-recipe",
    is_flag=True,
    help="Skip recipe execution even if --config or --recipe-path provides one.",
)
def cli(  # noqa: PLR0913
    config_path: Path | None,
    db_path: Path | None,
    recipe_path: Path | None,
    keys: tuple[str, ...],
    sample_size: int,
    expected_keys: tuple[str, ...],
    targets: tuple[str, ...],
    input_items: tuple[str, ...],
    db_field: str,
    skip_recipe: bool,  # noqa: FBT001
) -> None:
    """Validate LMDB structure and optionally run a recipe on sampled entries."""
    config_dict = load_config(config_path) if config_path is not None else {}
    resolved_db_path = _resolve_db_path(config_dict, db_path)
    resolved_recipe_path = _resolve_recipe_path(config_dict, recipe_path)
    resolved_expected_keys = _resolve_expected_keys(config_dict, expected_keys)
    extra_inputs = _parse_assignments(input_items)

    entry_count = _get_entry_count(resolved_db_path)
    sample_keys = list(keys) if keys else _get_sample_keys(resolved_db_path, sample_size)
    if not sample_keys:
        msg = f"No keys found in LMDB: {resolved_db_path}"
        raise click.ClickException(msg)

    click.echo(f"db_path={resolved_db_path}")
    click.echo(f"entry_count={entry_count}")
    click.echo(f"sample_keys={sample_keys}")
    if resolved_expected_keys is not None:
        click.echo(f"expected_top_keys={resolved_expected_keys}")

    failures = 0
    failures += _run_structure_check(
        db_path=resolved_db_path,
        keys=sample_keys,
        expected_top_keys=resolved_expected_keys,
    )

    if not skip_recipe and resolved_recipe_path is not None:
        failures += _run_recipe_check(
            sample_keys,
            RecipeTestConfig(
                db_path=resolved_db_path,
                recipe_path=resolved_recipe_path,
                db_field=db_field,
                targets=list(targets) if targets else None,
                config_dict=config_dict,
                extra_inputs=extra_inputs,
            ),
        )
    elif resolved_recipe_path is None:
        click.echo("\n[Recipe Test]\nskipped: no recipe configured")
    else:
        click.echo("\n[Recipe Test]\nskipped: --skip-recipe")

    if failures > 0:
        click.echo(f"\nResult: FAIL ({failures} issue(s))")
        sys.exit(1)

    click.echo("\nResult: PASS")


if __name__ == "__main__":
    cli()
