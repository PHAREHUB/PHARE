---
id: knowledge-in-repo
name: Non-obvious technique is explained
mode: warning
report: coderabbit
paths:
  - "src/**"
---

Look for non-obvious techniques the pull request introduces: an index shift or
offset, a ghost-cell width, an ordering constraint between steps, a numerical
trick, a workaround for SAMRAI or MPI behaviour, a magic constant.

Pass if each one is explained where a future reader will look: a comment at that
line, or the pull request description, or the documentation under doc/.
Fail if one is not explained anywhere. Name the line and say in one sentence what
the comment should state.
