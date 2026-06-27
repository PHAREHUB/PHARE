#!/usr/bin/env bash

# use like
#  mpirun -n $N ./tools/mpirun.sh python3 ...
# or
#  mpirun -n $N ./tools/mpirun.sh gdb -batch -ex run -ex bt --args python3 ...
# this will let you get a stack trace per rank
#

CWD="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" && cd $CWD/.. && CWD=$PWD

set -ex
RANK="$OMPI_COMM_WORLD_RANK"
PORT=$((5678 + RANK))
mkdir -p .log

(
  # do LD_PRELOAD HERE if required for ASAN
  # gdbserver localhost:$PORT $@
  perf record -o "perf.${RANK}.data" -Tga -F 1000 $@

  # uftrace -e record -d trace/_$RANK build/src/phare/phare-exe tests/functional/harris/harris_3d.py

) 1> >(tee $CWD/.log/${RANK}.out ) 2> >(tee $CWD/.log/${RANK}.err >&2 )

