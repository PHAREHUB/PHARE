#!/usr/bin/env bash
CWD="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" && cd $CWD
set -o pipefail
shopt -s expand_aliases
alias cls="clear; printf '\033[3J'"
cls

set -e

FILE="tests/core/numerics/ohm/test_tileset_ohm.cpp"
FILE="tests/core/numerics/ampere/test_tileset_ampere.cpp"
FILE="tests/core/data/tiles/test_tile.cpp"
FILE="tests/amr/data/field/test_fields_schedules.cpp"
FILE="tests/amr/level_initializer/test_hybrid_level_initializer.cpp"
FILE="tests/core/data/particles/test_particles_construction.cpp"
FILE="tests/core/numerics/faraday/test_tileset_faraday.cpp"
FILE="tests/core/numerics/ion_updater/test_multi_updater.cpp"
CARGS=${CARGS:-""}

set -x
(
  [ ! -f "bin/core/libphare_core.a" ] && (
    set -e
    mkn clean build -Kqp test_core -d google.test,+ -Oa "-fPIC" ${CARGS}
    mkn clean build -Kqp core -a "-fPIC" -x res/mkn/mpi -Og 0
  )

  [ ! -f "bin/amr/libphare_amr.a" ] && (
    set -e
    mkn clean build -Kqdp test_amr -a "-fPIC" -x res/mkn/mpi -Og 0
  )

  [[ $FILE == tests/core/* ]] && (
    set -e
    mkn clean build -p test_core -M "${FILE}" ${CARGS} -w mkn.gpu -x res/mkn/clang_cuda dbg "$@"
  ) || true

  [[ $FILE != tests/core/* ]] && (
    set -e
    mkn clean build -qp test_diagnostics -M "${FILE}" ${CARGS} -w mkn.gpu -x res/mkn/clang_cuda_mpi run "$@"
    echo "ok"
  ) || true

) 1> >(tee $CWD/.mkn.sh.out ) 2> >(tee $CWD/.mkn.sh.err >&2 )

# stop if not CI
[ -z "$CI" ] && exit 0 # continue for more

(
  set -xe
  mkn -tp core_tests ${CARGS} clean build test run "$@"
  mkn -tp more_tests ${CARGS} clean build test run "$@"
)
