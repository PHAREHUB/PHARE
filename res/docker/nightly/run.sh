#!/usr/bin/env bash
set -xe

eval "$(modulecmd bash load mpi/openmpi-x86_64)"

(
  cd /root
  git clone https://github.com/PHAREHUB/PHARE --depth 10 --recursive --shallow-submodules phare
  cd phare && mkdir build && cd build

  # we are not keeping these binaries, only output data - portability is not an issue
  CMAKE_CXX_FLAGS="-DNDEBUG -g0 -O3 -march=native -mtune=native"
  CMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS} -DPHARE_DIAG_DOUBLES=1"
  CMAKE_CONFIG="-DCMAKE_BUILD_TYPE=Release"
  cmake -G Ninja ${CMAKE_CONFIG} -DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS}" ..
  ninja
)

(
  cd /root/phare
  export PYTHONPATH="$PWD:$PWD/build:$PWD/pyphare"
  ./tools/data/generate.sh
)
rm -rf /root/phare
