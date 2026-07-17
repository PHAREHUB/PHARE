#!/usr/bin/env bash
CWD="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" && cd $CWD

set -e
. $CWD/mkn_func.sh

export PHARE_GPU_BYTES=500000000 # 500MB
export PHARE_SCOPE_TIMING=1
export PHARE_ASYNC_THREADS=1
export PHARE_COMPARE=0
# ASAN_OPTIONS=detect_leaks=0 OMPI_MCA_memory=^patcher

TEST="-M tests/core/numerics/ion_updater/test_multi_updater.cpp"
ARGS="${TEST} $(mkn_get_opts)"

(
  cd ..
  python3 -O tools/mkn_profile.py $ARGS $@
)
