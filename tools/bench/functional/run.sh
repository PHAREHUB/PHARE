#!/usr/bin/env bash
set -ex
CWD="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" && cd $CWD
ROOT=$(cd ../../.. && pwd)
N_CORES=${N_CORES:-1}

for F in $(find test_cases -maxdepth 1 -type f); do
  python3 $F;
done

(
    cd "$CWD"
    mkdir -p results
    for D in $(find "generated" -maxdepth 1 -type d | tail -n +2); do
        DIR="$( basename $D)"
        for F in $(find $D -maxdepth 1 -type f); do
            FILE=$(realpath -s --relative-to=$ROOT $F)
            echo "./build/src/phare/phare-exe" "$FILE" >>  "$CWD/results/${DIR}.sh"
        done
    done
)

(
    cd "$ROOT"
    python3 -m phlop.run.perf -d $CWD/results -i "*.sh" -c $N_CORES
)
