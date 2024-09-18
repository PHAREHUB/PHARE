#!/usr/bin/env bash
set -ex
CWD="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

(
    cd $CWD
    mpirun -n 4 tests/simulator/test_init_from_samrai.py
)
