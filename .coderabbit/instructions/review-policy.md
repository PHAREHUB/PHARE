---
id: review-policy
name: What to comment on in PHARE
---

PHARE is a hybrid PIC / MHD code with adaptive mesh refinement on SAMRAI. Its
reviewers come from physics, numerics, HPC and C++; write every comment so any
of them understands it without opening other files.

Comment on:
- wrong results: physics, numerics, indexing, centering, ghost cells, time levels;
- results that change with the number of MPI ranks, the patch layout or a restart;
- a helper, type or loop that already exists elsewhere in PHARE (name both places);
- a test that would still pass if the change were reverted;
- runtime cost added to code executed every time step.

Do not comment on, because the CI gate already checks them:
- formatting (clang-format), ruff findings, clang-tidy findings;
- `src/core` including SAMRAI, amr, simulator or python3 headers;
- removed tests, loosened test tolerances, missing PR body sections.

Format of each comment:
1. what the code is and does, in one plain line;
2. the problem as a concrete scenario: input or state, then what goes wrong;
3. how you know: read at file:line, or UNVERIFIED;
4. the fix, in one or two lines.
Mention a repeated pattern once, not on every occurrence.
