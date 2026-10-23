#!/usr/bin/env bash
CWD="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

CARGS=${CARGS:-""}
set -e

(
    cd $CWD/../../../../../..

    [ ! -f "bin/core/libphare_core.a" ] && (
        mkn clean build -tKOqp core -a "-fPIC" ${CARGS} -x res/mkn/mpi
        F=$(( ${#PWD} + 1 ))
        for f in $(find $CWD -name "*.cpp"); do
            A=${f:F}
            O=$(basename ${f%.*})
            mkn build -p test -tK -a -std=c++20 -M ${A} -w mkn.kul -o $O -O
        done
    )

    for f in $(find $CWD -name "*.cpp"); do
        A=${f:F}
        O=$(basename ${f%.*})
        mkn run -p test -M ${A} -o $O
    done
)
