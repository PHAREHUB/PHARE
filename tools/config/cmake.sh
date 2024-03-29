#!/usr/bin/env bash

set -e

CWD="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" && cd "$CWD"
ROOT_DIR=$(cd ../.. && pwd)
[ -f "${ROOT_DIR}/res/cmake/def.cmake" ] || (echo "Corruption detected" && exit 1) # confirm correct directory

[ -d build ] || (
    mkdir -p build && cd build
    cmake -DCMAKE_BUILD_TYPE=Release .. # 2>&1 > /dev/null
    make && ./phare_configurator
)
python3 config.py

# here as a placeholder while we flesh this out, but this does hookin to the normal cmake with -Dphare_configurator=ON
# probably to be generated by python (config.py) eventually
cat > local.cmake <<- EOM
cmake_minimum_required (VERSION 3.20.1)
project(configured_phare)
message("")
message("!!PHARE CONFIGURATED!!")
message("")
EOM
