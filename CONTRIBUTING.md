# Contributing to PHARE

## Issues

Every idea for changing the code — a bug fix or a new feature — starts as an
issue, before any code is written.

- **Bug reports** should make clear whether the bug is reproducible, and if
  so, how (steps, config, environment). If it isn't reproducible, describe
  what was observed and any conditions that seemed related.
- **Feature / change requests** should lead with the *why*: what problem this
  solves or what capability is missing. The *how* is secondary and can be
  worked out later, either in the issue itself or in the implementing PR.

Use the [bug report](.github/ISSUE_TEMPLATE/bug_report.yml) or
[feature request](.github/ISSUE_TEMPLATE/feature_request.yml) template when
opening an issue.

## Pull Requests

- **A PR should reference an issue.** PRs without an issue must be small,
  self-evident, and urgent — straightforward enough that why/priority/timing
  don't need discussion.
- **One issue per PR.** A PR should exist for a single reason to change the
  code. Bundling unrelated changes together slows down review — small,
  focused PRs merge faster.
- **Keep PRs small.** PRs of more than ~500 lines or more than ~10 files
  should be rare; split the work if possible.
- **PR description** should include:
  - A title that states the take-home message — why this PR, not just what
    it touches.
  - A link to the issue it implements, and an explanation of *how*: is the
    issue fully or partially implemented? What assumptions were made?
  - If the PR touches several files or is long, an explanation of the
    *design*. Diagrams (e.g. Mermaid sequence/class diagrams) are welcome
    for new classes or concepts.
- **Tests.** Bring unit or functional tests with the code as much as
  possible. If tests aren't included, explain why in the PR description.

Use the [PR template](.github/pull_request_template.md) — it mirrors this
checklist.
