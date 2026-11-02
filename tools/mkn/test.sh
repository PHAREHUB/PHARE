#!/usr/bin/env bash

set -e

CWD="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $CWD && cd ../.. && ROOT=$PWD

# export MKN_DBG="gdb -batch -ex run -ex bt --args"
# export MKN_DBG="/opt/rocm/bin/rocprofv2 --kernel-trace"
# export MKN_DBG="/opt/rocm/bin/rocprof -d outputFolder --hip-trace"
# export MKN_DBG="/opt/rocm/bin/rocgdb -batch -ex run -ex bt --args"

EMIT="-emit-llvm -S"
OPT_RECORD="" #-fsave-optimization-record=yaml -foptimization-record-file=f.yaml"
# then: /usr/lib/llvm-14/share/opt-viewer/opt-viewer.py f.yaml opts
OTHER="" #-DTHRUST_DEBUG_SYNC -ftemplate-backtrace-limit=0"
#CXXFLAGS="${CXXFLAGS} ${OPT_RECORD} ${OTHER}"

CXXFLAGS="-fPIC" # -ggdb -g -O0 -ftemplate-backtrace-limit=0"
XFIL="-x res/mkn/clang_cuda.yaml"
( which hipcc || [ -f /opt/rocm/bin/hipcc ] ) && XFIL="-x res/mkn/hip.yaml"

M_FILE="-M tests/core/data/particles/sorting/test_gpu_sorting_mkn.cpp"

# time (
  # export AMD_LOG_LEVEL=1
  # export AMD_SERIALIZE_KERNEL=3
  # export AMD_SERIALIZE_COPY=3
  # export HIP_TRACE_API=1
  # export HIP_PRINT_ENV=1

  # date
  
  mkn build ${M_FILE} -p test_core -a "${CXXFLAGS}" -w mkn.gpu ${XFIL} -Og 0 run
#   date
# ) 1> >(tee $ROOT/.mkn.sh.out ) 2> >(tee $ROOT/.mkn.sh.err >&2 )
