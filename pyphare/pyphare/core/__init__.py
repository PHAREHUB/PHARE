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


import dataclasses


@dataclasses.dataclass
class CliArgs:
    dry_run: bool = dataclasses.field(default_factory=lambda: False)
    write_reports: bool = dataclasses.field(default_factory=lambda: False)


def parse_cli_args():
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
        default=False,
    )

    return CliArgs(**vars(parser.parse_args()))
