#!/usr/bin/env bash

# use like
#  mpirun -n $N ./tools/mpirun.sh python3 ...
# or
#  mpirun -n $N ./tools/mpirun.sh gdb -batch -ex run -ex bt --args python3 ...
# this will let you get a stack trace per rank
#


CWD="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" && cd $CWD/.. && CWD=$PWD

RANK="$OMPI_COMM_WORLD_RANK"

mkdir -p .log

(
  # do LD_PRELOAD HERE if required for ASAN
  $@

) 1> >(tee $CWD/.log/${RANK}.out ) 2> >(tee $CWD/.log/${RANK}.err >&2 )

