#!/usr/bin/env bash

# use like
#  mpirun -n $N ./tools/mpirun.sh python3 ...
# or
#  mpirun -n $N ./tools/mpirun.sh gdb -batch -ex run -ex bt --args python3 ...
# this will let you get a stack trace per rank
#
#  You might need `OMPI_MCA_memory=^patcher mpirun -n ...` with ASAN
#


CWD="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" && cd $CWD/.. && CWD=$PWD

RANK="$OMPI_COMM_WORLD_RANK"

mkdir -p .log

(
  # do LD_PRELOAD HERE if required for ASAN
  # $MKN_DBG
  export PHARE_LOG=RANK_FILES
  export PHARE_SCOPE_TIMING=1
  export ASAN_OPTIONS=detect_leaks=0
  export LD_PRELOAD=/usr/lib64/libasan.so.8.0.0
  # $MKN_DBG
  $@

) 1> >(tee $CWD/.log/${RANK}.out ) 2> >(tee $CWD/.log/${RANK}.err >&2 )

