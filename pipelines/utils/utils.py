from importlib import import_module
from pathlib import Path
from typing import Any

from omegaconf import OmegaConf


def dotted_to_obj(path: str) -> object:
    """Convert a dotted path string to a Python object."""
    module_name, attr_name = path.rsplit(".", 1)
    module = import_module(module_name)
    return getattr(module, attr_name)


def load_config(config_path: Path) -> dict[str, Any]:
    """Load configuration from a YAML file using OmegaConf."""
    cfg = OmegaConf.load(config_path)
    cfg_dict = OmegaConf.to_container(cfg, resolve=True)
    if not isinstance(cfg_dict, dict):
        msg = "Configuration file must contain a dictionary at the top level."
        raise TypeError(msg)

    config: dict[str, Any] = {str(key): value for key, value in cfg_dict.items()}

    if "load_func" in config and isinstance(config["load_func"], str):
        config["load_func"] = dotted_to_obj(config["load_func"])
    if "transform_func" in config and isinstance(config["transform_func"], str):
        config["transform_func"] = dotted_to_obj(config["transform_func"])
    if "convert_func" in config and isinstance(config["convert_func"], str):
        config["convert_func"] = dotted_to_obj(config["convert_func"])
    if "project_func" in config and isinstance(config["project_func"], str):
        config["project_func"] = dotted_to_obj(config["project_func"])

    return config


def load_data_list(data_dir: Path, pattern: str = "*.cif*") -> list[Path]:
    """Load a list of data file paths from a directory."""
    return list(data_dir.rglob(pattern))
