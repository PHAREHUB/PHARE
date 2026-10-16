#!/usr/bin/env bash


set -e
CWD="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $CWD && cd ../.. && ROOT=$PWD # move to project root
ARGS=""
EXEC="build dbg"
DBOPT="-gO 0 "
X_FILE="-x mpi" #-x clang_asan"
FILES=(
  tests/amr/resources_manager/test_resources_manager.cpp
)
time (
  date
  for FILE in ${FILES[@]}; do
    mkn $EXEC -M ${FILE} -p test_amr -a "${ARGS}" $DBOPT $X_FILE $@
  done
  date
) 1> >(tee $ROOT/.mkn.sh.out ) 2> >(tee $ROOT/.mkn.sh.err >&2 )

###
#
# -r --gtest_filter=IonUpdaterTest/1.particlesUntouchedInMomentOnlyMode
###
