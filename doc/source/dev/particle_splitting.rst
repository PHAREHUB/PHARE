.. _dev-particle-splitting:

==================
Particle splitting
==================

.. container:: breadcrumb

   You are here: :doc:`/dev/index` › :doc:`/dev/one_time_step` ›
   :doc:`/dev/level_ghost_particles` › Particle splitting

**In one sentence:** to move plasma from a coarse :term:`level` to a finer one, PHARE
replaces each coarse particle by a few smaller particles placed around it, with
the same velocity, so that they deposit almost the same density on the fine grid.

Why splitting
-------------

In PHARE, a particle's shape (the small cloud it deposits on the grid) is
measured in cells of the grid it lives on. A coarse particle moved to a fine level,
even with its weight rescaled to carry the same plasma (see the weight trick
below), would see its cloud shrink to fine-cell size, half as wide in each
direction. The density would become lumpy, and the
fine level would have :math:`2^d` times fewer particles per cell than the coarse
one (:math:`d` the dimension). Splitting it into :math:`N` children, spread
around the parent's position, rebuilds almost the same wide cloud out of
fine-grid clouds, and gives the fine level more particles per cell.

Splitting is used every time fine particles must be made from coarse ones:

* to fill a **new or rebuilt refined level** (its domain particles);
* to fill the **level ghost layer** around a refined level: Old when a level
  is created or rebuilt, New at the first fine substep of each coarse step
  (:doc:`/dev/level_ghost_particles`).

The picture
-----------

One coarse particle (thick line) and its children (dashed) on the fine grid,
in 1D. The grey curve is the sum of the children: it is close to the parent.
Distances are in fine cells.

.. plot::
   :include-source: False
   :caption: The density deposited by one coarse particle (black) and by its
      split children (grey: their sum). Order 1 and 2 are the particle shapes
      used for interpolation orders 1 and 2. Values from split_1d.hpp.

   import numpy as np
   import matplotlib.pyplot as plt

   def B(k, u):
       u = np.abs(u)
       if k == 1:
           return np.clip(1 - u, 0, None)
       if k == 2:
           return np.where(u < 0.5, 0.75 - u**2,
                           np.where(u < 1.5, 0.5 * (1.5 - u)**2, 0.0))
       return np.where(u < 1, 2/3 - u**2 + u**3 / 2,
                       np.where(u < 2, (2 - u)**3 / 6, 0.0))

   # (order, label, [(weight, delta), ...]) from split_1d.hpp
   cases = [
       (1, "order 1, N = 2", [(0.5, 0.551569), (0.5, -0.551569)]),
       (1, "order 1, N = 3", [(0.5, 0.0), (0.25, 1.0), (0.25, -1.0)]),
       (2, "order 2, N = 2", [(0.5, 0.663959), (0.5, -0.663959)]),
   ]
   x = np.linspace(-4, 4, 801)
   fig, axes = plt.subplots(3, 1, figsize=(6, 8), sharex=True)
   for ax, (k, title, children) in zip(axes, cases):
       ax.plot(x, B(k, x / 2) / 2, color="k", lw=2.5, label="coarse particle")
       total = np.zeros_like(x)
       for w, d in children:
           c = w * B(k, x - d)
           total += c
           ax.plot(x, c, ls="--", color="#2a7ab9", lw=1.2)
           ax.plot([d], [0], marker="^", color="#2a7ab9", ms=8, clip_on=False)
       ax.plot(x, total, color="0.55", lw=1.5, label="sum of children")
       ax.set_title(title, fontsize=12)
       ax.set_xticks(range(-4, 5))
       ax.tick_params(labelsize=11)
       ax.grid(axis="x", color="0.9")
       for s in ("top", "right"):
           ax.spines[s].set_visible(False)
   axes[0].legend(fontsize=11, frameon=False, loc="upper left")
   axes[-1].set_xlabel("fine cells", fontsize=12)
   fig.tight_layout()

Where the children go in 2D and 3D with the default number of children
(blue: children, black: parent). In 2D they sit on the diagonals, in 3D on the
axes.

.. plot::
   :include-source: False
   :caption: Default split patterns, interpolation order 1. 2D: 4 children at
      (±0.572, ±0.572) fine cells, weight 1/4 each. 3D: 6 children at ±0.966
      fine cells along each axis, weight 1/6 each (drawn in the x-y plane, the
      two along z are on top of the parent). Values from split_2d.hpp and
      split_3d.hpp.

   import matplotlib.pyplot as plt

   fig, axes = plt.subplots(1, 2, figsize=(6, 3.4))
   d2, d3 = 0.571783, 0.966431
   pts2 = [(-d2, -d2), (-d2, d2), (d2, -d2), (d2, d2)]
   pts3 = [(d3, 0), (-d3, 0), (0, d3), (0, -d3)]
   for ax, pts, title in ((axes[0], pts2, "2D, N = 4"), (axes[1], pts3, "3D, N = 6")):
       for v in range(-2, 3):
           ax.axvline(v, color="0.85", lw=0.8)
           ax.axhline(v, color="0.85", lw=0.8)
       ax.plot([0], [0], "o", color="k", ms=10)
       ax.plot([p[0] for p in pts], [p[1] for p in pts], "o", color="#2a7ab9", ms=8)
       ax.set_xlim(-2, 2); ax.set_ylim(-2, 2); ax.set_aspect("equal")
       ax.set_title(title, fontsize=12)
       ax.set_xticks(range(-2, 3)); ax.set_yticks(range(-2, 3))
       ax.tick_params(labelsize=11)
       ax.set_xlabel("fine cells", fontsize=12)
   fig.tight_layout()

The flow inside the refine operator:

.. mermaid::

   flowchart TD
     S["coarse domain particles<br/>(local, or streamed from another rank)"] --> F["move to fine index space<br/>toFineGrid"]
     F --> Q{"inside the destination box<br/>grown by splitBoxGrowth?"}
     Q -- no --> X["skipped"]
     Q -- yes --> P["split into N children<br/>Splitter / PatternDispatcher"]
     P --> K{"child inside the<br/>destination box?"}
     K -- no --> Y["dropped"]
     K -- yes --> A["append to the destination array:<br/>domain, or level ghost Old / New"]
     click A href "level_ghost_particles.html" "Level ghost particles"

How it works
------------

The work is done by ``ParticlesRefineOperator``, the :term:`SAMRAI` refine operator
PHARE provides for particles
(:src:`src/amr/data/particles/refine/particles_data_split.hpp#L54-L281`).
SAMRAI calls its ``refine`` method with a coarse source patch, a fine destination
patch and the list of destination boxes to fill.

#. **Source particles.** Only the coarse *domain* particles are split (L166).
   When the coarse :term:`patch` is on another :term:`MPI` :term:`rank`, SAMRAI
   first streams the coarse
   particles of the region of interest into a temporary patch, and they arrive
   in its domain array (comment L159-L165). The region SAMRAI gathers around
   the destination is set by ``getStencilWidth``: the particle ghost width, in
   coarse cells (L80-L84).
#. **To fine index space.** ``toFineGrid`` rewrites the particle position
   (cell index ``iCell`` plus offset ``delta`` in [0, 1)) in fine cells:
   ``iCell * 2 + floor(delta * 2)``, and the remainder as the new ``delta``
   (L33-L48).
#. **Candidates.** A particle is split only if this fine position is inside the
   destination box grown by ``splitBoxGrowth`` fine cells (``getSplitBox``,
   L267-L280, test at L215). Parents just outside the box can still put a child
   inside it.
#. **Split.** The ``Splitter`` makes the :math:`N` children
   (``PatternDispatcher::dispatch``,
   :src:`src/amr/data/particles/refine/splitter.hpp#L85-L119`). Each child gets
   the parent's velocity, the parent's position plus an offset ``δ`` (in fine
   cells), and a weight (next section).
#. **Keep.** Children inside the destination box are appended to the array
   chosen at compile time by the operator's split type (L225-L259):

   .. list-table::
      :header-rows: 1

      * - Split type
        - Destination array
        - Used when
      * - ``interior``
        - ``domainParticles``
        - a refined level is created or rebuilt
      * - ``coarseBoundaryOld``
        - ``levelGhostParticlesOld``
        - level created or rebuilt
      * - ``coarseBoundaryNew``
        - ``levelGhostParticlesNew``
        - first fine substep of each coarse step

   The enum (L25-L30) has a fourth value, ``coarseBoundary``, which writes into
   ``levelGhostParticles``; no refiner uses it at this commit (the three
   operators actually built are listed in ``RefinementParams``, L287-L302).

The tricks
----------

Weights are multiplied by :math:`2^d`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** A particle weight in PHARE is a density per particle, not an
amount of plasma: at initialization, ``weight = n / nppc`` (density over
particles per cell, ``cellWeight`` in
:src:`src/core/data/ions/particle_initializers/maxwellian_particle_initializer.hpp#L219`).
A fine cell has :math:`2^d` times less volume than a coarse cell (:math:`d` the
dimension). The same amount of plasma in a smaller cell is a larger density.

**Fix.** A child's weight is the parent weight × its pattern weight ×
:math:`2^d` (``fineParticle.weight``,
:src:`src/amr/data/particles/refine/splitter.hpp#L104`). The pattern weights
of each splitter sum to 1, so the children carry exactly the parent's plasma.

The split box is grown by the largest child offset
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** Children are placed up to ``δ`` fine cells away from their parent.
If the code only looked at parents inside the destination box, children of
parents just outside it would be missing, and the density at the box edge would
be too low.

**Fix.** ``splitBoxGrowth`` grows the box by the larger of half the particle
size and the largest offset rounded up (``PatternDispatcher::maxDelta``)
(:src:`src/amr/data/particles/refine/splitter.hpp#L182-L190`). It matters for
the 2D, order 1, 8-children pattern, whose diagonal children reach 1.058 fine
cells.

**History.** #1369 (splitters checked against reference tables).

Patterns and weights come from published tables
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** The children's offsets and weights decide how well the fine density
matches the coarse one. They are not obvious to derive.

**Fix.** Each ``Splitter`` specialization combines a few *patterns* (one child
in the centre, children on the axes, on the diagonals...), each with one weight
and one offset (``SplitPattern``,
:src:`src/amr/data/particles/refine/splitter.hpp#L35-L45`). The values in
``split_1d.hpp``, ``split_2d.hpp`` and ``split_3d.hpp`` follow the tables of
`arXiv:2104.10675 <https://arxiv.org/abs/2104.10675>`_ (Computer Physics
Communications 261, 107666, 2021). A unit test checks the weights and offsets
against those tables.

**History.** #1369 (splitters checked against reference tables); before it,
some non-default splitters did not conserve the plasma (#1340).

The number of children is fixed at compile time
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The number of children :math:`N` is a template parameter: each
(dimension, interpolation order, :math:`N`) combination is a separately
compiled simulator. Users choose :math:`N` with ``refined_particle_nbr`` in the
simulation block. The Python modules built by default are those listed in
:src:`res/sim/all.txt#L1-L33` (read by :src:`src/python3/CMakeLists.txt#L5-L9`;
a build can point ``PHARE_PERMUTATIONS`` to another list). The default
:math:`N` (bold) is the first value of each row of ``valid_refined_particle_nbr``
(:src:`pyphare/pyphare/pharein/simulation.py#L206-L215`, used at L236-L237).

.. list-table:: Combinations built by default (``res/sim/all.txt``)
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

.. note::

   In 3D, the Python check ``valid_refined_particle_nbr`` accepts more values
   than the default build provides (for example 8 or 64); with those, no
   compiled simulator exists.

Go deeper
---------

This is the bottom of this branch. The unit tests for the splitters are in
``tests/amr/data/particles/refine/``.

Related
-------

* :doc:`/dev/level_ghost_particles` — the main user of splitting during a run.
* :doc:`/dev/field_refinement` — the field counterpart: refinement operators for
  E and B.
* :doc:`/numerics/particles` — particle shapes and interpolation orders.
