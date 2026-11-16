#!/usr/bin/env bash

set -e

CWD="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" && cd "$CWD"
ROOT_DIR=$(cd ../.. && pwd)
[ -f "${ROOT_DIR}/res/cmake/def.cmake" ] || (echo "Corruption detected" && exit 1) # confirm correct directory

export CMAKE_COMMAND=${1}
export CMAKE_CXX_COMPILER=${2}
export PYTHON_EXECUTABLE=${3}
export PHARE_BUILD_DIR=${4}

PROBE_DIR="${PHARE_BUILD_DIR}/configurator"
PROBE_OUTPUT="${PROBE_DIR}/PHARE_MPI_Get_library_version.txt"

if [ ! -f "${PROBE_OUTPUT}" ]; then
    rm -rf "${PROBE_DIR}"
    mkdir -p "${PROBE_DIR}"
    if ! (
        cd "${PROBE_DIR}" &&
        "$CMAKE_COMMAND" -DCMAKE_BUILD_TYPE=Release "${CWD}" \
                         -DCMAKE_CXX_COMPILER="${CMAKE_CXX_COMPILER}" &&
        "$CMAKE_COMMAND" --build . --verbose &&
        ./phare_configurator &&
        [ -f "${PROBE_OUTPUT}" ]
    ); then
        rm -rf "${PROBE_DIR}"
        echo "PHARE configurator: the probe in ${PROBE_DIR} failed to build or run" >&2
        exit 1
    fi
fi

"$PYTHON_EXECUTABLE" config.py || { rm -rf "${PROBE_DIR}"; exit 1; }
