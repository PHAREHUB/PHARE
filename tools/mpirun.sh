#!/usr/bin/env bash

# use like
#  mpirun -n $N ./tools/mpidbg.sh python3 script.py
#
# available options
#  gdbserver: for parallel debugging like with vscode
#  perf: for performance metrics sampling
#  gdb: for printing stacktraces on any rank (default)

# ./tools/mpirun.sh python3 -m phlop.run.test_cases -i tests/simulator/ -r Initialization2DTest.test_levelghostparticles_have_correct_split_from_coarser_particle_14 --prefix="$MKN_DBG"

CWD="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

set -ex
RANK="${OMPI_COMM_WORLD_RANK:-0}"
PORT=$((5678 + RANK))

cd "$CWD/.." && CWD=$PWD


[ ! -f ".asan.options" ] && cat >.asan.options <<EOL
interceptor_via_fun:*__cxa_throw*
interceptor_via_fun:*__interception::real___cxa_throw*

EOL

(
  mkdir -p .log
  ## for asan ignore weird errors in python/ft2font
  LIB_ASAN=$(readlink -e $(g++ -print-file-name=libasan.so))
  LIB_STD=$(readlink -e $(g++ -print-file-name=libstdc++.so))
  ASAN_OPTIONS=suppressions=.asan.options:detect_leaks=0 LD_PRELOAD=$LIB_ASAN:$LIB_STD "$@"

  ## for debugging mpirun in parallel processes
  # gdbserver localhost:$PORT "$@"

  ## for sampling function times
  # perf record -o "perf.${RANK}.data" -Tga -F 1000 "$@"

  ## stop/print on first error
  # gdb -batch -ex run -ex bt --args "$@"

) 1> >(tee "$CWD/.log/${RANK}.out" ) 2> >(tee "$CWD/.log/${RANK}.err" >&2 )

