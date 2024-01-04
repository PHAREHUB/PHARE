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

    return parser.parse_args()
