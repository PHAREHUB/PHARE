#!/usr/bin/env bash

# use like
#  mpirun -n $N ./tools/mpidbg.sh python3 script.py
#
# available options
#  gdbserver: for parallel debugging like with vscode
#  perf: for performance metrics sampling
#  gdb: for printing stacktraces on any rank (default)

CWD="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

set -ex
RANK="$OMPI_COMM_WORLD_RANK"
PORT=$((5678 + RANK))

(
  cd "$CWD/.." && CWD=$PWD
  mkdir -p .log

  # do LD_PRELOAD HERE if required for ASAN

  ## for debugging mpirun in parallel processes
  # gdbserver localhost:$PORT "$@"

  ## for sampling function times
  # perf record -o "perf.${RANK}.data" -Tga -F 1000 "$@"

  ## stop/print on first error
  gdb -batch -ex run -ex bt --args "$@"

) 1> >(tee "$CWD/.log/${RANK}.out" ) 2> >(tee "$CWD/.log/${RANK}.err" >&2 )

