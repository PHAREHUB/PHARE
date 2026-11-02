#!/usr/bin/env bash
CWD="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

CARGS=${CARGS:-""}
set -e

(
    cd $CWD/../../../../../..

    for f in $(find $CWD -name "*.cpp"); do
        A=${f:F}
        O=$(basename ${f%.*})
        mkn run -p test -M ${A} -o $O
    done

) 1> >(tee $CWD/.run.sh.out ) 2> >(tee $CWD/.run.sh.err >&2 )

