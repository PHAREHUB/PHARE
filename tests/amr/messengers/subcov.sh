#!/usr/bin/env bash
#
# Author: Philip Deegan
# Email : philip.deegan@gmail.com
# Date  : 28 - October - 2019
#
# Usage:
#  Scans gmon.out and greps for strings, string not found = error
#
######################################################################

set -ex

CWD="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
IFS=$'\n'
TEST="test-messenger"


TO_GREP_LIST=(
  PHARE::amr::FieldData<PHARE::core::GridLayout.*::unpackStream
)

for TO_GREP in "${TO_GREP_LIST[@]}"; do
  WC=$(grep "TO_GREP". gprof.* | wc -l)
  $(( WC == 0 )) && \
    echo "Coverage tests failed:" && \
    echo "String: $TO_GREP" && \
    echo "Does not exist for test $TEST" && \
    exit 1
done

