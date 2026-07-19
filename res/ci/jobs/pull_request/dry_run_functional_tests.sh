#!/usr/bin/env bash

(
  set -ex

  # prevent openmpi complaining about root user on CI
  export OMPI_ALLOW_RUN_AS_ROOT=1
  export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

  cd build
  cmake .. -DPHARE_EXEC_LEVEL_MIN=11 -DPHARE_EXEC_LEVEL_MAX=100
  export PHARE_DRY_RUN=1
  ctest --output-on-failure
)
