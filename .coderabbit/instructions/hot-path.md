---
id: hot-path
name: Runtime cost on the per-step path
paths:
  - "src/core/numerics/**"
  - "src/core/data/particles/**"
  - "src/amr/messengers/**"
  - "src/amr/solvers/**"
  - "src/amr/multiphysics_integrator.hpp"
---

This code runs every time step, often per particle or per cell. Flag:
- heap allocation, `std::function`, or virtual dispatch inside a per-particle or
  per-cell loop;
- copies of fields or particle arrays where a view or reference works;
- new MPI collectives or barriers, and whether they run once or every step;
- new template instantiations or heavy includes in widely included headers.
Give an estimate with its arithmetic (for example, per particle per step) and
label it "estimate, not measured". Do not present an estimate as a measurement.
