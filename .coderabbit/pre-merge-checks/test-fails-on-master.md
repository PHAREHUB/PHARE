---
id: test-fails-on-master
name: Behaviour change comes with a test
mode: warning
report: coderabbit
paths:
  - "src/**"
  - "pyphare/pyphare/**"
---

Decide whether the pull request changes behaviour: results, outputs, accepted
inputs, or a fixed bug. Pure refactors, renames, comments and performance changes
that claim identical results do not change behaviour.

Pass if behaviour does not change; or if a test added or changed in this pull
request exercises the new behaviour and would fail on the target branch without
the source change (name the test and the assertion); or if the description's
Tests section explains why no test is possible.
Fail otherwise. Name the behaviour that is untested and suggest the smallest test
that would fail on the target branch.
