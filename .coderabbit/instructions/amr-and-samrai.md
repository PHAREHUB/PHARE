---
id: amr-and-samrai
name: AMR, SAMRAI and parallel correctness
paths:
  - "src/amr/**"
  - "src/simulator/**"
---

Check:
- several quantities in one SAMRAI refine algorithm: quantities of the same
  patch-data type but different centering must not fall into one equivalence
  class (SAMRAI compares types, not centerings, and would compute wrong overlaps);
- schedules: created once per regrid, not inside the time loop;
- ghost cells: the ghost width covers the stencil reach; fills at level borders
  use the right time interpolation;
- periodic boundaries: shifts handled, patches touching only at a corner included;
- MPI: every rank reaches every collective; no result depends on message arrival
  order, rank count or patch decomposition; nothing writes to shared state from
  one rank only;
- regrid: data moved to new patches, not recomputed from stale values;
- particle splitting (`src/amr/data/particles/refine/`): children copy the
  parent's velocity; fine weight = parent weight x pattern weight x refinement
  ratio^dim (`splitter.hpp`), and the pattern weights sum to 1 to float precision
  (stored as `float`, e.g. `split_1d.hpp`), so charge and momentum are conserved
  to about 1e-6 relative, not to double round-off; children are kept only where
  the operator claims (fine domain or level ghost); the result does not depend on
  patch layout or particle order.
