.. _dev-particle-splitting:

==================
Particle splitting
==================

.. container:: breadcrumb

   You are here: :doc:`/dev/index` › :doc:`/dev/one_time_step` ›
   :doc:`/dev/level_ghost_particles` › Particle splitting

**In one sentence:** to move plasma from a coarse :term:`level` to a finer one,
PHARE replaces each coarse particle by N smaller particles placed around it,
with the same velocity, so that they deposit almost the same density on the
fine grid.

Why: a particle's shape is measured in cells of its own grid. Moved unchanged
to the fine grid, its cloud would shrink to half its width and the density would
get lumpy. N children spread around it rebuild nearly the same wide cloud.

The children rebuild the parent's cloud
---------------------------------------

.. plot::
   :include-source: False
   :caption: Density deposited by one coarse particle (black) and by its
      children (dashed; grey: their sum), in fine cells. Order 1 with N = 3 is
      exact. Values from split_1d.hpp.

   import numpy as np
   from phare_draw import *
   def S(k, u):
       u = np.abs(u)
       if k == 1:
           return np.clip(1 - u, 0, None)
       return np.where(u < 0.5, 0.75 - u**2, np.where(u < 1.5, 0.5 * (1.5 - u)**2, 0.0))
   cases = [(1, "Order 1, N = 2", [(0.5, 0.551569), (0.5, -0.551569)]),
            (1, "Order 1, N = 3", [(0.5, 0.0), (0.25, 1.0), (0.25, -1.0)]),
            (2, "Order 2, N = 2", [(0.5, 0.663959), (0.5, -0.663959)])]
   fig, axes = strip(3, [c[1] for c in cases], height=1.6, aspect=False)
   x = np.linspace(-4, 4, 801)
   for a, (k, _, kids) in zip(axes, cases):
       a.plot(x, 4 * S(k, x / 2) / 2, color=C["ink"], lw=2.5)
       tot = np.zeros_like(x)
       for w, d in kids:
           c = 4 * w * S(k, x - d); tot += c
           a.plot(x, c, ls="--", color=C["L1"], lw=1.2)
           a.plot([d], [0], marker="^", color=C["L1"], ms=8)
       a.plot(x, tot, color=C["grey"], lw=1.5)
       for v in range(-4, 5):
           text(a, v, -0.35, str(v), color=C["grey"])
   finish(fig, axes)

* Weight of a child = parent weight × its pattern weight × 2^d; the children's
  weights sum to the parent's × 2^d (pattern weights sum to 1, to float
  precision) (:src:`src/amr/data/particles/refine/splitter.hpp#L99-L117`).
* 2D default: 4 children on the diagonals at ±0.572; 3D default: 6 on the axes
  at ±0.966 fine cells (order 1; ``split_2d.hpp``, ``split_3d.hpp``).

Inside the refine operator
--------------------------

.. plot::
   :include-source: False
   :caption: Making fine particles for a destination box (blue) from one coarse
      particle, order 1, N = 4. Thin grid: fine cells; thick: coarse cells.
      δ: position inside the cell, from 0 to 1 in each direction.

   from phare_draw import *
   from matplotlib.patches import Rectangle
   dest = (6, 0, 12, 6)
   p = (5.6, 2.8)
   d = 0.571783
   kids = [(p[0] + sx * d, p[1] + sy * d) for sx in (-1, 1) for sy in (-1, 1)]
   def base(a, grown=False):
       grid(a, (2, -1, 13, 7), color=C["light"])
       for v in range(2, 14, 2):
           a.plot([v, v], [-1, 7], color=C["grey"], lw=1.2, zorder=1)
       for v in range(0, 7, 2):
           a.plot([2, 13], [v, v], color=C["grey"], lw=1.2, zorder=1)
       a.add_patch(Rectangle(dest[:2], 6, 6, fc=C["patch"], ec=C["patch_edge"], lw=2, zorder=0.5))
       if grown:
           a.add_patch(Rectangle((5, -1), 8, 8, fc="none", ec=C["patch_edge"], lw=1.5, ls="--", zorder=2))
   fig, (a1, a2, a3) = strip(3, [
       "Coarse particle in fine index space (toFineGrid)",
       "Inside the box grown by splitBoxGrowth: split",
       "Keep only the children inside the destination box"], height=2.4)
   base(a1); particles(a1, [p], color=C["ink"], size=60)
   text(a1, 3.6, 5.6, "coarse cell (2,1)\nδ = (0.8, 0.4)", color=C["ink"])
   base(a2, grown=True); particles(a2, [p], color=C["ink"], size=60)
   text(a2, 11.0, 7.4, "grown by 1 fine cell", color=C["patch_edge"])
   particles(a2, kids, color=C["L1"], size=45)
   base(a3); particles(a3, [p], color=C["ink"], size=60, alpha=0.25)
   particles(a3, [k for k in kids if k[0] >= 6], color=C["ok"], size=45)
   for k in kids:
       if k[0] < 6:
           text(a3, k[0], k[1], "✕", color=C["bad"], bold=True)
   text(a3, 9, 3, "kept", color=C["ok"], bold=True); text(a3, 3.4, 3, "dropped", color=C["bad"], bold=True)
   finish(fig, [a1, a2, a3])

* Only coarse *domain* particles are split; on another rank SAMRAI first streams
  them into a temporary patch (``ParticlesRefineOperator::refine_``,
  :src:`src/amr/data/particles/refine/particles_data_split.hpp#L155-L166`).
* ``toFineGrid``: ``iCell × 2 + floor(delta × 2)``, remainder as new delta
  (:src:`src/amr/data/particles/refine/particles_data_split.hpp#L33-L48`).
* Split if the fine position is in the destination box grown by
  ``splitBoxGrowth`` (L203-L215, ``getSplitBox`` L267-L280).
* Children inside the destination box go to the array of the operator's type:
  ``interior`` → domain (new or rebuilt level), ``coarseBoundaryOld`` /
  ``coarseBoundaryNew`` → level ghost Old / New (L225-L259; ``RefinementParams``
  L287-L302; a fourth type, ``coarseBoundary``, is unused).

The tricks
----------

Weights are multiplied by 2^d
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** A weight is a density per particle (``n / nppc`` at init,
:src:`src/core/data/ions/particle_initializers/maxwellian_particle_initializer.hpp#L219`);
a fine cell has 2^d times less volume.
**Fix.** ``fineParticle.weight = weight × pattern weight × 2^d``
(:src:`src/amr/data/particles/refine/splitter.hpp#L104`).

The split box is grown by the largest child offset
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** Parents just outside the box put children inside it (2D order 1
N = 8 reaches 1.058 fine cells).
**Fix.** ``splitBoxGrowth`` = max(half the particle size, ceil(max δ))
(:src:`src/amr/data/particles/refine/splitter.hpp#L182-L190`). :pr:`1369`.

Patterns come from published tables
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Offsets and weights follow `arXiv:2104.10675 <https://arxiv.org/abs/2104.10675>`_
(CPC 261, 107666, 2021), checked by a unit test. :pr:`1369` (splitters checked
against tables); before it, some non-default splitters lost plasma (:issue:`1340`).

N is fixed at compile time
^^^^^^^^^^^^^^^^^^^^^^^^^^

``refined_particle_nbr`` picks a compiled simulator; the default build compiles
:src:`res/sim/all.txt#L1-L33` (``PHARE_PERMUTATIONS`` can point elsewhere,
:src:`src/python3/CMakeLists.txt#L5-L9`). Default = first value of each row
(:src:`pyphare/pyphare/pharein/simulation.py#L206-L215`, taken at L237).

.. list-table:: N built by default (bold: default)
   :header-rows: 1

   * - Dimension
     - Order 1
     - Order 2
     - Order 3
   * - 1D
     - **2**, 3
     - **2**, 3, 4
     - **2**, 3, 4, 5
   * - 2D
     - **4**, 5, 8, 9
     - **4**, 5, 8, 9, 16
     - **4**, 5, 8, 9, 25
   * - 3D
     - **6**, 12, 27
     - **6**, 12, 27
     - **6**, 12, 27

N = 2^d keeps the coarse number of particles per cell (defaults in 1D, 2D); the
3D default gives 6/8 of it. In 3D the Python check accepts values that are not
built (e.g. 8, 64).

Go deeper
---------

Bottom of this branch. Unit tests: ``tests/amr/data/particles/refine/``.

Related
-------

* :doc:`/dev/level_ghost_particles` · :doc:`/dev/field_refinement` ·
  :doc:`/numerics/particles`
