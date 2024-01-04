"""
PHARE build and runtime validation checks
"""

import os
import sys
import json
import dataclasses
from pathlib import Path

from pyphare.core import phare_utilities
from pyphare import cpp

DOT_PHARE_DIR = Path(os.getcwd()) / ".phare"


ERROR_MESSAGES = dict(
    on_phare_config_error="""
        Warning: PHARE was not built with the configurator active""",
    on_python3_version_mismatch="""
        Inconsistency detected!
        Python during build and now are not the same!
        Build python version  : {}
        Current python version: {}""",
    on_build_config_check_runtime_error="""
        Could not interrogate python versions
        Please see 'Having Issues' Section of ISSUES.TXT
        Actual error (consider adding to issue report): {}""",
)


def print_error(key, *args):
    print(ERROR_MESSAGES[key].format(*args).strip())


def python_version_from(binary):
    return phare_utilities.decode_bytes(
        phare_utilities.run_cli_cmd(f"{binary} -V", check=True).stdout.strip()
    )


def check_build_config_is_runtime_compatible(strict=True):
    try:
        build_config: dict = cpp.build_config()

        if "PHARE_CONFIG_ERROR" in build_config:
            print_error("on_phare_config_error")
            return

        build_python_version = build_config["PYTHON_VERSION"]
        current_python_version = python_version_from(sys.executable)
        if build_python_version != current_python_version:
            print_error(
                "on_python3_version_mismatch",
                build_python_version,
                current_python_version,
            )

            raise ValueError("Python version mismatch!")

    except RuntimeError as e:
        print_error("on_build_config_check_runtime_error", e)

    except ValueError as e:
        print(e)
        if strict:
            raise e


@dataclasses.dataclass
class RuntimeSettings:
    python_version: str
    python_binary: str


def try_system_binary(cli, log_to):
    with open(log_to, "w") as f:
        try:
            proc = phare_utilities.run_cli_cmd(cli, check=True)
            f.write(phare_utilities.decode_bytes(proc.stdout).strip())
        except Exception as e:
            f.write(f"failed to run cli command {cli}\n")
            f.write(f"error {e}")


def try_system_binaries(log_dir):
    try_system_binary("free -g", log_dir / "free_dash_g.txt")
    try_system_binary("lscpu", log_dir / "lscpu.txt")
    try_system_binary("hwloc-info", log_dir / "hwloc_info.txt")


def log_runtime_config():
    cpp_lib = cpp.cpp_lib()

    settings = RuntimeSettings(
        python_binary=sys.executable, python_version=python_version_from(sys.executable)
    )

    if cpp_lib.mpi_rank() == 0:
        DOT_PHARE_DIR.mkdir(exist_ok=True, parents=True)
    cpp_lib.mpi_barrier()

    rank_dir = DOT_PHARE_DIR / f"rank_{cpp_lib.mpi_rank()}"
    rank_dir.mkdir(exist_ok=True)

    with open(rank_dir / "runtime_config.json", "w") as f:
        json.dump(dataclasses.asdict(settings), f)

    try_system_binaries(rank_dir)
