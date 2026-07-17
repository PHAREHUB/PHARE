#!/usr/bin/env bash
CWD="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

CARGS=${CARGS:-""}
set -e

(
    cd $CWD/../../../../../..

    mkn clean build -tKOqp core -a "-fPIC" ${CARGS} -x res/mkn/mpi -g 0
    F=$(( ${#PWD} + 1 ))
    i=0
    pids=()
    for f in $(find $CWD -name "*.cpp"); do
        A=${f:F}
        O=$(basename ${f%.*})
        mkn build -p test -tKa "-std=c++20 -DPHARE_FORCE_LOG_LINE=1" -M ${A} -o $O -Ow mkn.avx -g 0 &
        pids[${i}]=$!
        i+=1
    done

    for pid in ${pids[*]}; do
        wait $pid
    done

) 1> >(tee $CWD/.build.sh.out ) 2> >(tee $CWD/.build.sh.err >&2 )
