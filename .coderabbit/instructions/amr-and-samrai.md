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
- regrid: data moved to new patches, not recomputed from stale values.
