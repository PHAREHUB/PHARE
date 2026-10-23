#!/usr/bin/env bash
set -exu
CWD="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" && cd $CWD && cd ../.. && ROOT=$PWD

shopt -s expand_aliases
alias cls="clear; printf '\033[3J'"

MKN_ARGS="-P mkn.base=gpu_ -Ox res/mkn/hip_mpi " # -W

run(){
    FILE="$1"
    mkn clean build -M "${FILE}" run -p test_core ${MKN_ARGS}
}

run_mpi(){
    FILE="$1"
    N="$2"
    MKN_DBG="mpirun -n 3" mkn clean build run -M "${FILE}" -p test_core ${MKN_ARGS}
}

. "$CWD/_run_gpu_tests.sh"

