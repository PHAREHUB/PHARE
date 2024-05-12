#!/usr/bin/env bash

set -e

CWD="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" && cd "$CWD"
ROOT_DIR=$(cd ../.. && pwd)
[ -f "${ROOT_DIR}/res/cmake/def.cmake" ] || (echo "Corruption detected" && exit 1) # confirm correct directory

export CMAKE_COMMAND=${1}
export CMAKE_CXX_COMPILER=${2}
export PYTHON_EXECUTABLE=${3}

[ -d build ] || (
    mkdir -p build && cd build
    "$CMAKE_COMMAND" -DCMAKE_BUILD_TYPE=Release .. \
                     -DCMAKE_CXX_COMPILER="${CMAKE_CXX_COMPILER}" # 2>&1 > /dev/null
    make VERBOSE=1 && ./phare_configurator
)
$PYTHON_EXECUTABLE config.py
