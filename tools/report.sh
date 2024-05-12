#!/usr/bin/env bash

# usage:
#  ./tools/report.sh
#    generates zip archive for logging issues

set -eu
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
(
    export PYTHONPATH="${PWD}:${PWD}/build:${PWD}/pyphare:${PYTHONPATH}"
    python3 "$ROOT/tools/report.py"
)
