---
id: numerics-and-physics
name: Numerics and physics
paths:
  - "src/core/numerics/**"
  - "src/core/data/**"
  - "src/amr/solvers/**"
  - "src/amr/data/**"
  - "pyphare/pyphare/pharein/**"
---

Check, for changed formulas and stencils:
- centering: each quantity is read and written at its Yee-grid location
  (primal or dual per direction, from `GridLayout`); a stencil mixing centerings
  without interpolation is a bug;
- time level: the predictor-corrector reads fields at the time level it needs
  (n, n+1/2, n+1); flag a read of a quantity that was already advanced;
- conservation: mass, momentum, energy and div B are kept where the scheme claims
  it; name the term that breaks it;
- accuracy order: a change of stencil or interpolation states the order before
  and after;
- normalisation: PHARE works in normalised units; flag a hard-coded physical
  constant or a mixed unit;
- refinement and coarsening: a quantity refined from coarse to fine and coarsened
  back stays consistent (sums for moments, averages for fields as the operator
  defines).
When a doubt can be settled by a small numerical test (convergence order,
conservation over a few steps), propose that test.
