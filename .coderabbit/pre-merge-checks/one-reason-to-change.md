---
id: one-reason-to-change
name: One reason to change
mode: warning
report: coderabbit
---

List the distinct logical changes in the pull request: feature, bug fix,
refactor, formatting pass, unrelated cleanup, dependency bump.

Two changes are independent when each could be reverted or merged alone without
breaking the other. A large change that serves one purpose, such as a rename
touching many files, counts as one change.

Pass if the pull request contains one change, or several that depend on each
other.
Fail if it contains two or more independent changes. Name each one, list the
files or hunks that belong to it, and say how to split the pull request.
