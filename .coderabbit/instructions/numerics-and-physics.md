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
  (primal or dual per direction, from `GridLayout`). A derivative flips the
  centering of its direction, so terms built from differently centered fields
  can be combined (see `CurlComponent` in `differential_operators.hpp`). The bug
  is two terms that land on different centerings and are added without
  interpolating one of them;
- time level: the predictor-corrector reads fields at the time level it needs
  (n, n+1/2, n+1); flag a read of a quantity that was already advanced;
- conservation: mass, momentum, energy and div B are kept where the scheme claims
  it; name the term that breaks it;
- accuracy order: a change of stencil or interpolation states the order before
  and after;
- normalisation: PHARE works in normalised units; flag a hard-coded physical
  constant or a mixed unit;
- MHD (finite volume / Godunov): wave-speed estimates (`mhd_speeds.hpp`) are not
  smaller than the true fastest speed used by Rusanov, HLL and HLLD; a
  reconstruction keeps its claimed order, and a TVD limiter (minmod, van Leer)
  does not create new extrema (WENO and MP5 may overshoot slightly); constrained
  transport keeps div B at round-off (flag an average that breaks it); primitive
  and conservative conversions are inverse to round-off; density and pressure
  stay positive (say where a floor is applied);
- refinement and coarsening: a quantity refined from coarse to fine and coarsened
  back stays consistent (sums for moments, averages for fields as the operator
  defines).
When a doubt can be settled by a small numerical test (convergence order,
conservation over a few steps), propose that test.
