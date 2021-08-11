#! /usr/bin/env bash

set -ex

CWD="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

mkdir -p $CWD/build
[ ! -d  "$CWD/build" ] && echo "mkn expects cmake configured build directory" && exit 1

MKN_X_FILE=${MKN_X_FILE:-settings}
MKN_MPI_X_FILE=${MKN_MPI_X_FILE:-res/mkn/mpi}
MKN_GPU_X_FILE=${MKN_GPU_X_FILE:-res/mkn/clang_nvcc}
export MKN_LIB_LINK_LIB=1

exit 0

# verify compiler setup
# mkn clean build -x ${MKN_GPU_X_FILE} -ga -DKUL_GPU_CUDA dbg -p mkn.gpu_depositor run test

[ ! -d ./subprojects/thrust ] && \
  git clone https://github.com/NVIDIA/thrust \
    -b main subprojects/thrust --depth 10 --shallow-submodules --recursive
[ ! -d ./subprojects/cuda-samples ] && \
  git clone https://github.com/NVIDIA/cuda-samples \
    -b main subprojects/cuda-samples --depth 10 --shallow-submodules --recursive

 
# gtest doens't like mpi compilers
#mkn clean build -x ${MKN_X_FILE} -p test_diagnostics -tKOd google.test,+

#mkn clean build -x ${MKN_MPI_X_FILE} -dtKOp py

# mkn clean build -x ${MKN_MPI_X_FILE} -p test_diagnostics -KO 9 test run
# mkn clean build -x ${MKN_GPU_X_FILE} -Oa -DKUL_GPU_CUDA run

#mkn clean build -x ${MKN_GPU_X_FILE} -Op cpp
