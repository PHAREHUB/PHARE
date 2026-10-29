#!/usr/bin/env bash
CWD="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$CWD"/..

set -e

PY_FILES=$(find . -name "*.py")

python3 -m black phlop tests
pylint --errors-only phlop tests
isort phlop tests

for FILE in ${PY_FILES[@]}; do

  autoflake -i "$FILE"


done
