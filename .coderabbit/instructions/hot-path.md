---
id: hot-path
name: Runtime and memory cost
paths:
  - "src/core/numerics/**"
  - "src/core/data/particles/**"
  - "src/amr/data/**"
  - "src/amr/messengers/**"
  - "src/amr/solvers/**"
  - "src/amr/load_balancing/**"
  - "src/amr/multiphysics_integrator.hpp"
  - "src/diagnostic/**"
  - "pyphare/pyphare/pharesee/**"
---

Most of this code runs every time step, often per particle or per cell; the rest
runs at each regrid, dump or read of a large run. Flag:
- heap allocation, `std::function`, or virtual dispatch inside a per-particle or
  per-cell loop;
- copies of fields or particle arrays where a view or reference works;
- new MPI collectives or barriers, and whether they run once or every step;
- new template instantiations or heavy includes in widely included headers;
- memory: new per-cell or per-particle storage, and objects kept alive longer
  than before (a lifetime change is easy to miss);
- code outside the time step that still costs on large runs: regrid and load
  balancing, diagnostics writing, and pharesee reading large runs; say when it
  runs (every step, every regrid, every dump, once) and what it scales with;
- a claimed speed-up or an "identical cost" claim without a measurement: ask for
  the case, the build and the before/after numbers;
- new timing or profiling code: check first what SAMRAI timers
  (`tbox::TimerManager`, the `xfer::` and `tbox::Schedule::` timers) and phlop
  already provide, and ask for per-rank imbalance, not only totals.
Give an estimate with its arithmetic (for example, per particle per step) and
label it "estimate, not measured". Do not present an estimate as a measurement.
