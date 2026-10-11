#!/usr/bin/env bash
# ruff (with fixes) and black over all python the CI gate checks, see ruff.toml
#  for staged lines only, see res/git_hooks/pre-commit
set -euo pipefail

cd "$(git rev-parse --show-toplevel)"

paths=(pyphare tests tools)
rc=0
ruff check "${paths[@]}" --fix "$@" || rc=$? # black runs either way
python3 -m black "${paths[@]}"
exit $rc
