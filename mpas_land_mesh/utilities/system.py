import sys
import os
from pathlib import Path
from typing import Tuple


def get_python_environment() -> Tuple[str, str, str]:
    """
    Retrieve the Python environment path, name, and type from the interpreter.

    Supports conda, venv, virtualenv, pyenv, and system Python installations.

    Returns
    -------
    Tuple[str, str, str]
        A tuple containing:
        - env_path : str
            Absolute path to the Python environment directory.
        - env_name : str
            Name of the environment (directory name, or 'system' for system Python).
        - env_type : str
            Type of environment: 'conda', 'venv', 'virtualenv', 'pyenv', or 'system'.

    Raises
    ------
    RuntimeError
        If unable to determine Python executable path, or if the detected
        environment path does not exist or is not a directory.
    """
    python_path = sys.executable

    if not python_path:
        raise RuntimeError("Unable to determine Python executable path")

    python_path_obj = Path(python_path).resolve()

    # Check for conda environment
    # Method 1: Check for 'envs' directory in path
    components = str(python_path_obj).split(os.sep)
    if "envs" in components:
        try:
            envs_index = components.index("envs")
            if envs_index + 1 < len(components):
                env_name = components[envs_index + 1]
                env_path = os.sep.join(components[: envs_index + 2])
                env_type = "conda"

                env_path_obj = Path(env_path)
                if not env_path_obj.exists() or not env_path_obj.is_dir():
                    raise RuntimeError(
                        f"Detected conda environment path '{env_path}' does not exist or is not a directory"
                    )

                return env_path, env_name, env_type
        except (ValueError, IndexError):
            pass  # Fall through to other checks

    # Method 2: Check CONDA_DEFAULT_ENV variable (for conda base or named envs)
    conda_env = os.environ.get("CONDA_DEFAULT_ENV")
    if conda_env:
        env_path = sys.prefix
        env_name = conda_env
        env_type = "conda"

        env_path_obj = Path(env_path)
        if env_path_obj.exists() and env_path_obj.is_dir():
            return env_path, env_name, env_type

    # Check for venv (has pyvenv.cfg file)
    current = python_path_obj.parent
    for _ in range(3):  # Check up to 3 levels up
        pyvenv_cfg = current / "pyvenv.cfg"
        if pyvenv_cfg.exists():
            env_path = str(current)
            env_name = current.name
            env_type = "venv"
            return env_path, env_name, env_type
        current = current.parent

    # Check for virtualenv (VIRTUAL_ENV variable)
    virtual_env = os.environ.get("VIRTUAL_ENV")
    if virtual_env:
        env_path = virtual_env
        env_name = Path(virtual_env).name
        env_type = "virtualenv"

        env_path_obj = Path(env_path)
        if env_path_obj.exists() and env_path_obj.is_dir():
            return env_path, env_name, env_type

    # Check for pyenv
    # Method 1: Check for 'versions' directory in path
    if "versions" in components:
        try:
            versions_index = components.index("versions")
            if versions_index + 1 < len(components):
                env_name = components[versions_index + 1]
                env_path = os.sep.join(components[: versions_index + 2])
                env_type = "pyenv"

                env_path_obj = Path(env_path)
                if env_path_obj.exists() and env_path_obj.is_dir():
                    return env_path, env_name, env_type
        except (ValueError, IndexError):
            pass

    # Method 2: Check PYENV_VERSION variable
    pyenv_version = os.environ.get("PYENV_VERSION")
    if pyenv_version:
        env_path = sys.prefix
        env_name = pyenv_version
        env_type = "pyenv"

        env_path_obj = Path(env_path)
        if env_path_obj.exists() and env_path_obj.is_dir():
            return env_path, env_name, env_type

    # Default to system Python
    env_path = str(python_path_obj.parent)
    env_name = "system"
    env_type = "system"

    env_path_obj = Path(env_path)
    if not env_path_obj.exists() or not env_path_obj.is_dir():
        raise RuntimeError(
            f"Python executable directory '{env_path}' does not exist or is not a directory"
        )

    return env_path, env_name, env_type


def get_extension_from_path(filepath: str) -> str:
    """
    Extract the file extension from a given filename.

    Args:
        filename (str): The name of the file.

    Returns:
        str: The file extension (including the dot), or an empty string if no extension is found.
    """
    basename = os.path.basename(filepath)
    _, ext = os.path.splitext(basename)
    return ext