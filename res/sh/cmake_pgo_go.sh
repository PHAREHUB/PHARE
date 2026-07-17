#!/usr/bin/env bash
# usage:
#  build PHARE with profile guilded optimizations

set -ex
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $SCRIPT_DIR && cd .. && CWD=$PWD # move to project root

[ ! -f "$CWD/CMakeLists.txt" ] && echo "script expected to be run from project root" && exit 1
[ ! -d "$CWD/subprojects/cppdict/include" ] && git submodule update --init

export OPYTHONPATH=${PYTHONPATH}
export PYTHONPATH=${PWD}:${PWD}/build:${PWD}/pyphare:${OPYTHONPATH}

THREADS=${THREADS:-"10"}
BUILD_DIR=${BUILD_DIR:-"$CWD/build"}
SAMRAI=${SAMRAI:-""} # "" = from system or as subproject if not found in system
[[ -n "${SAMRAI}" ]] && SAMRAI=-DSAMRAI_ROOT="${SAMRAI}"

CMAKE_CXX_FLAGS="-DNDEBUG -g0 -O3 -march=native -mtune=native"
CMAKE_CONFIG=" -Dphare_configurator=ON -DCMAKE_BUILD_TYPE=Release"

export CC=${CC:-"gcc"}
export CXX=${CXX:-"g++"}
set -xe

mkdir -p ${BUILD_DIR}
(
    cd ${BUILD_DIR}

    cmake $CWD ${SAMRAI} ${CMAKE_CONFIG} -G Ninja -DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS}" -DPGO_GEN=ON;
    ninja -v -j${THREADS} cpp cpp_etc dictator
    python3 tests/functional/harris/harris_2d_2.py

    cmake $CWD ${SAMRAI} ${CMAKE_CONFIG} -G Ninja -DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS}" -DPGO_GEN=OFF -DPGO_USE=ON;
    ninja -v -j${THREADS} cpp cpp_etc dictator
    python3 tests/functional/harris/harris_2d_2.py # should be faster than without any PGO
)
