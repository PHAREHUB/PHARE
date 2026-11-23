---
id: expressiveness
name: Readable C++ inside functions
paths:
  - "src/**/*.hpp"
  - "src/**/*.cpp"
---

Apply in this order, to changed functions only:
1. One level of abstraction: state what the function name promises, list the
   steps that explain it; each statement maps to one step. A line that only
   makes an earlier step work (index arithmetic, a counter, a buffer resize) is a
   sub-step: suggest a named helper. Length alone is not a finding.
2. Names: a reader guesses what a function, variable or condition does from its
   name. Flag `data`, `info`, `manager`, `temp`, bare `true`/`false` arguments at
   call sites, and conditions that need a comment (suggest a named bool).
3. Illegal states: flag same-typed parameters that are easy to swap, sentinel
   values (`-1`, `nullptr` meaning absent) where `std::optional` or a strong type
   (like `DimConst<>`) fits.
4. Control flow: long if / else-if chains that hide a rule; suggest named
   conditions or early returns.
5. Raw loops doing what a standard algorithm or a PHARE utility names.
6. A function whose name hides that it mutates state.
Label these comments "expressiveness" and keep them below correctness findings.
At most 3 expressiveness comments per review: the ones a reader would gain most
from.
