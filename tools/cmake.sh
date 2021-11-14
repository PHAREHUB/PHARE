#!/usr/bin/env bash

# To ignore this file from git : git update-index --skip-worktree tools/cmake.sh
#  to re-enable git indexing (changes) for this file use: git update-index --no-skip-worktree tools/cmake.sh
#

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $SCRIPT_DIR && cd .. && CWD=$PWD # move to project root
exec 19>$CWD/.cmake.sh.cmd # set -x redirect
export BASH_XTRACEFD=19  # set -x redirect
THREADS=${THREADS:="$(nproc --all)"}
BUILD_DIR=${BUILD_DIR:="$CWD/build"}
SAMRAI=${SAMRAI:=""} # "" = as subproject
FFF=("${BUILD_DIR}")
set -xe
time (
  date
  [ -n "$CLEAN" ] && (( $CLEAN == 1 )) && for f in ${FFF[@]}; do rm -rf $f; done
  [ ! -f "$CWD/CMakeLists.txt" ] && echo "script expected to be run from project root" && exit 1
  [ ! -d "$CWD/subprojects/cppdict/include" ] && git submodule update --init
  mkdir -p ${BUILD_DIR}
  [[ -n "${SAMRAI}" ]] && SAMRAI=-DSAMRAI_ROOT="${SAMRAI}"
  (cd ${BUILD_DIR} && cmake $CWD ${SAMRAI} -G Ninja -DdevMode=ON -Dasan=OFF \
    -DCMAKE_CXX_FLAGS="-g3 -O3 -march=native -mtune=native -DPHARE_DIAG_DOUBLES=1")
  ninja -C ${BUILD_DIR} -v -j${THREADS}
  date
) 1> >(tee $CWD/.cmake.sh.out ) 2> >(tee $CWD/.cmake.sh.err >&2 )
