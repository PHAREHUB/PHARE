# Red explainer

Explains red CI builds to the PR author: what failed, why, and a fix written as a prompt to paste
into a coding agent. Advisory only: it never blocks a merge.

1. `tc.py` reads TeamCity through guest access (no token).
2. `triage.py` sorts each red by fixed rules: the CI's fault (registry, timeout, false red) gets
   one line and no AI; the PR's fault (compile error, failed test) goes on.
3. `excerpt.py` cuts the few log lines that matter.
4. `agent.py` runs a free OpenCode model in a folder holding the excerpt, the PR diff and the
   built commit's code; it may only read (`opencode.json`: no edits, no shell, no web).
5. `comment.py` lays out one comment per PR, each answer inside a code block: GitHub renders
   no link, HTML or mention there, whatever the model wrote.

`explain.py` holds what replay and CI share (sorting, git, one agent run).

## In CI: `.github/workflows/ci-explainer.yml`

TeamCity reports each PR build to GitHub as a commit status. GitHub runs the workflow on every
status, always with master's workflow and scripts.

1. `ci.py check SHA`: nothing while a TeamCity status on the commit is pending; then `red`, or
   `green` if an earlier comment needs correcting. Skips old commits (the PR moved), commits
   already reported, and, during the trial, PRs without the `ci-explain` label.
2. `ci.py explain`: the agent, in a job with a read-only token. At most 4 distinct failures per
   commit. Its inputs and raw output stay 14 days as the workflow artifact `ci-explainer`.
3. `ci.py post` / `ci.py green`: the only job allowed to write; one comment per PR, edited, and only
   if the commit is still the PR's head (an explanation takes minutes).

Run by hand: Actions → ci-explainer → Run workflow, with the PR head's full sha.

The model is chosen once, `MODEL` in `agent.py`. To switch it without a PR, set the repository
variable `CI_EXPLAINER_MODEL` (Settings → Secrets and variables → Actions → Variables), e.g.
`opencode/<model>-free`; delete it to go back to the default. Free models may train on what they
read: they only ever see public code and logs.

## Replay (offline, nothing posted)

```
python3 res/ci/explainer/replay.py 20261019 /tmp/replay              # all reds since that day
python3 res/ci/explainer/replay.py 20261019 /tmp/replay --no-agent   # sorting and excerpts only
python3 res/ci/explainer/replay.py 20261019 /tmp/replay --only 100729,99827
```

Writes `would_post.md` (the comments), `grading.md` (one row per explained failure, to mark
right / useless / wrong) and `summary.json`. Needs `opencode` on the PATH and an `upstream` remote
for PHAREHUB/PHARE.

## Tests

```
python3 -m unittest discover -s res/ci/explainer/tests
```
