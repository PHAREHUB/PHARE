#!/usr/bin/env bash
#
# Author: Philip Deegan
# Email : philip.deegan@gmail.com
# Date  : 28 - October - 2019
#
# Usage:
#  Scans subdirectories for subcov.sh and gmon.out
#   Executes subcov.sh to validate call graph
#   Example subcov.sh is in tests/.../messengers
#   Must be executed from CMake build dir, or project root
#
#   subcov.sh files are executed for isolation
#    run: find tests -name subcov.sh | xargs chmod +x
#
######################################################################

set -ex

BUILD_=${BUILD:-$PWD}
CWD="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
pushd $CWD/.. 2>&1 > /dev/null
ROOT=$PWD
popd 2>&1 > /dev/null

# IFS = Internal Field Separato. See "globbing"
IFS=$'\n'

pushd $CWD 2>&1 > /dev/null
for f in $(find . -name subcov.sh); do
  [[ $f == "./"* ]] && f=${f:2}
  DIR=tests/$(dirname $f)
  [ -f $BUILD/$DIR/gmon.out ] \
    && pushd $BUILD/$DIR 2>&1 > /dev/null \
    && $CWD/$f && popd 2>&1 > /dev/null
done
popd 2>&1 > /dev/null
