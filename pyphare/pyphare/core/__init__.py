#
#
#
# program help looks like
"""
usage: phare_sim.py [-h] [-d]

options:
  -h, --help     show this help message and exit
  -d, --dry-run  Validate but do not run simulations
"""

import sys
import dataclasses


def disabled_for_testing():
    # check if any module is loaded from PHARE tests directory
    from pathlib import Path

    test_dir = Path(__file__).resolve().parent.parent.parent.parent / "tests"
    if test_dir.exists():
        test_dir = str(test_dir)
        for k, mod in sys.modules.items():
            if hasattr(mod, "__file__") and mod.__file__ and test_dir in mod.__file__:
                return True
    return False


@dataclasses.dataclass
class CliArgs:
    dry_run: bool = dataclasses.field(default_factory=lambda: False)
    write_reports: bool = dataclasses.field(default_factory=lambda: False)


def parse_cli_args():
    default_off = len(sys.argv) == 1 and disabled_for_testing()
    if default_off:
        return CliArgs()

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-d",
        "--dry-run",
        help="Validate but do not run simulations",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "-w",
        "--write_reports",
        help="Write build and runtime configs to disk",
        action="store_true",
        default=True,
    )

    return CliArgs(**vars(parser.parse_args()))
