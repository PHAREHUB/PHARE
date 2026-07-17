#!/usr/bin/env bash
set -ex

(
  cd /github/workspace
  mkdir build && cd build
  cmake -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON \
        -DlowResourceTests=ON            \
        -DCMAKE_CXX_FLAGS="-DPHARE_DIAG_DOUBLES=1 -O2" ..
  make VERBOSE=1
  ctest -j 2 --output-on-failure
)
