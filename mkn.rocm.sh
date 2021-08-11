#! /usr/bin/env bash

set -ex

CWD="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

MKN_X_FILE=${MKN_X_FILE:-settings}
MKN_MPI_X_FILE=${MKN_MPI_X_FILE:-res/mkn/mpi}
MKN_GPU_X_FILE=${MKN_GPU_X_FILE:-res/mkn/hcc}

mkdir -p $CWD/build
[ ! -d  "$CWD/build" ] && echo "mkn expects cmake configured build directory" && exit 1

export MKN_LIB_LINK_LIB=1

# verify compiler setup
# mkn clean build -x ${MKN_GPU_X_FILE} -ga -DKUL_GPU_ROCM dbg -p mkn.gpu_depositor run test

# git clone https://github.com/NVIDIA/thrust -b 1.13.0 subprojects/thrust --depth 10 --shallow-submodules --recursive

# gtest doens't like mpi compilers
mkn clean build -x ${MKN_X_FILE} -p test_diagnostics -tKOd google.test,+

mkn clean build -x ${MKN_MPI_X_FILE} -dtKOp py

mkn clean build -x ${MKN_MPI_X_FILE} -p test_diagnostics -KO 9 test run

mkn clean build -x ${MKN_GPU_X_FILE} -Oa -DKUL_GPU_ROCM run
