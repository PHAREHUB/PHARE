.. _dev-grabbing-particles:

=========================================================
Grabbing particles from neighbours; patch ghost particles
=========================================================

.. container:: breadcrumb

   You are here: :doc:`/dev/index` › :doc:`/dev/one_time_step` ›
   :doc:`/dev/ppc_order` › Grabbing particles from neighbours

**In one sentence:** a particle that leaves its patch for a neighbour patch of
the same level is first kept by the patch it left, as a *patch ghost particle*,
so that this patch can deposit it; then the neighbour *grabs* it into its own
domain.

Why it works this way
---------------------

A :term:`level` is cut into :term:`patches <patch>`, often on different
:term:`MPI` :term:`ranks <rank>`. Each patch pushes only its own particles. After a push,
some particles sit in a cell that belongs to a neighbour. Two things must
happen:

* the particle must **change owner**: from the next step on, the neighbour
  pushes it;
* its **density and flux** at the new time must still be counted, once.

PHARE does both in a fixed order: the patch the particle left deposits it, then
the particle is handed over. The neighbour does not deposit it in this step: the
nodes they share are completed later by a sum over patch borders
(:doc:`/dev/moments`).

The picture
-----------

.. plot::
   :include-source: False
   :caption: Two patches of the same level. A particle of patch A (black dot)
      is pushed into patch B. It lands in A's particle ghost box, in the part
      covered by B (blue): A keeps it as a patch ghost particle and deposits it.
      Then B grabs it: B takes from A's patch ghost array every particle in
      A's box grown by the particle ghost width, intersected with B's domain.
      Width 2 cells here (interpolation orders 2 and 3; 1 cell at order 1).

   import matplotlib.pyplot as plt
   from matplotlib.patches import Rectangle, FancyArrowPatch

   fig, ax = plt.subplots(figsize=(5.5, 3.8))
   g = 2
   A = (2, 2, 8, 8)
   B = (10, 2, 8, 8)
   for P, name in ((A, "patch A"), (B, "patch B")):
       x0, y0, w, h = P
       for x in range(x0, x0 + w + 1):
           ax.plot([x, x], [y0, y0 + h], color="0.8", lw=0.6)
       for y in range(y0, y0 + h + 1):
           ax.plot([x0, x0 + w], [y, y], color="0.8", lw=0.6)
       ax.add_patch(Rectangle(P[:2], w, h, fill=False, lw=2, color="k"))
       ax.text(x0 + w / 2, y0 + h + 0.4, name, ha="center", fontsize=14)
   ax.add_patch(Rectangle((B[0], B[1]), g, B[3], color="#7fb7e6", alpha=0.8))
   ax.add_patch(Rectangle((A[0] - g, A[1] - g), A[2] + 2 * g, A[3] + 2 * g,
                          fill=False, ls="--", lw=1.5, color="#333333"))
   ax.plot([8.6], [6.3], "o", color="k", ms=9, alpha=0.35)
   ax.plot([10.7], [6.8], "o", color="k", ms=9)
   ax.add_patch(FancyArrowPatch((8.6, 6.3), (10.6, 6.75), arrowstyle="->",
                                mutation_scale=15, color="k"))
   ax.text(12.3, 7.0, "after the push:\npatch ghost\nof A", fontsize=13, bbox=dict(fc="white", ec="none", alpha=0.9))
   ax.text(0, -1.4, "A's particle ghost box (dashed)", fontsize=13)
   ax.set_xlim(-0.8, 18.5); ax.set_ylim(-1.8, 12.6)
   ax.set_aspect("equal"); ax.axis("off")
   fig.tight_layout()

The order of operations in the second, real push of a time step
(the first push only predicts moments and exchanges nothing, see
:doc:`/dev/ppc_order`):

.. mermaid::

   sequenceDiagram
     participant A as patch A
     participant B as patch B
     A->>A: push domain particles
     A->>A: leavers in a neighbour's cells → A's patch ghost array
     A->>A: leavers in the level ghost layer → erased
     A->>A: deposit domain + patch ghost particles
     A-->>B: sum density and flux on shared nodes
     B->>A: grab: copy A's patch ghosts that are in B's domain
     A->>A: clear patch ghost array
     B->>B: drop particles in cells B does not own

How it works
------------

#. **Push and sort.** ``IonUpdater::updateAndDepositAll_`` pushes the domain
   particles. A *selector* is a test on the particle's cell used to sort the
   array; here ``inOwnedDomain`` puts at the end of the array the particles that
   left the cells owned by the patch
   (:src:`src/core/numerics/ion_updater/ion_updater.hpp#L258-L265`).
#. **Keep the leavers that went to a neighbour.** Among the leavers, those in
   the patch's *non-level-ghost box* (its domain plus the parts of its particle
   ghost box covered by same-level neighbours, i.e. everything except the level
   ghost layer;
   ``makeNonLevelGhostBoxFor``,
   :src:`src/amr/resources_manager/amr_utils.hpp#L224-L243`) are appended to the
   patch ghost array (L265-L271). All leavers are then removed from the domain
   array (L289). Leavers that went into the level ghost layer are lost on
   purpose: the coarse level represents that region
   (:doc:`/dev/level_ghost_particles`).
#. **Deposit.** Domain particles *and* patch ghost particles are deposited on
   the patch's grid (L306-L309). A leaver thus adds its share to the nodes it
   touches, including ghost nodes that this patch shares with the neighbour.
#. **Sum moments over patch borders.** The solver sums density and flux on the
   nodes shared by neighbour patches (``fillFluxBorders``,
   ``fillDensityBorders``, :src:`src/amr/solvers/solver_ppc.hpp#L620-L627`;
   :doc:`/dev/moments`).
#. **Grab.** ``HybridHybridMessengerStrategy::fillIonGhostParticles``
   (:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L375-L385`)
   runs a :term:`SAMRAI` :term:`schedule` over the level (``RefinerType::ExteriorGhostParticles``,
   :src:`src/amr/messengers/refiner.hpp#L196-L199`). For each pair of
   neighbours, the :term:`fill pattern` ``ParticleDomainFromGhostFillPattern`` defines
   the region to move: the source patch box grown by the particle ghost width,
   intersected with the destination domain
   (:src:`src/amr/data/particles/particles_variable_fill_pattern.hpp#L80-L90`).
   ``ParticlesData`` then copies the source's *patch ghost* particles in that
   region into the destination's *domain* array: ``copy_from_ghost`` on the same
   rank, ``pack_from_ghost`` / ``unpack_from_ghost`` across ranks
   (:src:`src/amr/data/particles/particles_data.hpp#L709-L793`).
#. **Clear.** Right after the schedule, every patch ghost array is emptied
   (L381-L383). Patch ghost particles live only between a push and the next
   exchange.
#. **Single owner.** Each patch drops the domain particles that landed in cells
   it does not own (``Ions::eraseDomainParticlesIn``,
   :src:`src/amr/solvers/solver_ppc.hpp#L651-L656`). This only matters when
   patches overlap (next section).

The tricks
----------

Deposit before handing over
^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** If the neighbour deposited the particle after receiving it, the
patch it left would also need its contribution on its own ghost and border
nodes, and the order between deposit, exchange and sum would get tangled.

**Fix.** The patch the particle left deposits it, from its patch ghost array,
in the same call as its domain particles
(:src:`src/core/numerics/ion_updater/ion_updater.hpp#L306-L309`). The border
sum then adds that contribution to the neighbour's copy of the shared nodes.
The particle is counted once, by one patch.

The exchange region is the source box grown by the ghost width
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** A patch ghost particle is at most a particle ghost width away from
the patch it left. The neighbour must take all of those, and only those, that
fall in its domain, also when patches overlap.

**Fix.** The region is ``grow(source box, particle_ghost_width) ∩ destination
domain`` (:src:`src/amr/data/particles/particles_variable_fill_pattern.hpp#L84-L90`).
A ``static_assert`` checks that SAMRAI's particle data ghost width covers that
distance, since SAMRAI only builds a transfer between patches that close
(L60-L63).

**History.** :pr:`1334` (overlap single owner).

Overlapping patches: one owner per cell
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** SAMRAI can produce patches of the same level that overlap. A
particle in the overlap would be held by both patches, pushed and deposited
twice.

**Fix.** Each cell has one owner: the patch with the smallest global id among
those that cover it (``makeNonOwnedBoxes``,
:src:`src/amr/resources_manager/amr_utils.hpp#L246-L280`). Each patch computes,
without communication, the boxes it does not own; the push treats them like a
neighbour's cells (``isOwned``,
:src:`src/core/numerics/ion_updater/ion_updater.hpp#L135-L143`), and after the
grab the particles that arrived there are dropped
(:src:`src/amr/solvers/solver_ppc.hpp#L651-L656`). Without overlap the
non-owned set is empty and nothing changes.

**History.** :pr:`1334` (overlap single owner).

Periodic neighbours
^^^^^^^^^^^^^^^^^^^

With periodic boundaries, a patch on one side of the domain is a neighbour of a
patch on the other side. SAMRAI gives the exchange a *transformation* (a shift
by the domain size); ``copy_from_ghost`` and ``pack_from_ghost`` shift the cell
index of each particle by that offset
(:src:`src/amr/data/particles/particles_data.hpp#L713-L730`).

How far a particle can move in one half push
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** All of the above assumes leavers stay close to the patch they
left.

**Fix.** The Boris pusher moves positions in two half steps. In each, it
computes the new offset from the start of the particle's current cell,
``delta + v Δt / (2 Δx)`` (``delta`` in [0, 1)), and stops the run with an error
if its size exceeds 2 (``MoveTwoCellException``,
:src:`src/core/numerics/pusher/boris.hpp#L157-L172`). So a particle moves at
most 2 − ``delta`` cells forward and 2 + ``delta`` cells backward (up to about
3) in one half step. When this happens the run aborts; the terminal only shows
``MPI_ABORT was invoked``, and the message ``PrePush`` (or ``PostPush``)
``Particle moved 2 cells`` is written in ``.log/<rank>.out``. The particle arrays
index particles by cell over a margin wide enough for any allowed move
(``particleMapMargin``, :src:`src/core/data/grid/gridlayout.hpp#L59-L68`).

**History.** :pr:`1592` (selectors missed far particles).

.. admonition:: Known issue
   :class: note

   The list of same-level neighbours of a patch is built with SAMRAI's default
   search width of 1 cell (``makeSameLevelNeighbors``,
   :src:`src/amr/resources_manager/amr_utils.hpp#L216-L220`), while the
   particle ghost width is 2 cells at interpolation orders 2 and 3. In 2D and
   3D, a neighbour that touches a patch's ghost box only near a corner can be
   missed; a particle pushed more than one cell into that corner is then
   treated as a level ghost and erased. A fix is proposed in
   `#1593 <https://github.com/PHAREHUB/PHARE/pull/1593>`_ (corner neighbour
   search).

Go deeper
---------

* :doc:`/dev/moments` — the sum over patch borders that completes the shared
  nodes.

Related
-------

* :doc:`/dev/level_ghost_particles` — the coarse-fine counterpart: particles
  entering a refined level from the coarser one.
* :doc:`/dev/ppc_order` — where the grab sits in a time step.
* :doc:`/dev/particles_ions` — the particle arrays of a patch.
* :doc:`/dev/ghost_cells` — patch ghosts and level ghosts for fields.
