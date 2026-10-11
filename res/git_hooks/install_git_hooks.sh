#!/usr/bin/env bash
# points git at the hooks in this directory, so hooks added here later need no reinstall
# undo: git config --unset core.hooksPath
set -euo pipefail

cd "$(git rev-parse --show-toplevel)"

hooks=res/git_hooks
current=$(git config --get core.hooksPath || true)

if [ -n "$current" ] && [ "$current" != "$hooks" ]; then
  echo "core.hooksPath is already set to '$current', not overriding" >&2
  exit 1
fi

git config core.hooksPath "$hooks"
echo "git hooks installed from $hooks"
