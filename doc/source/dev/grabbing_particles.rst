.. _dev-grabbing-particles:

=========================================================
Grabbing particles from neighbours; patch ghost particles
=========================================================

.. container:: breadcrumb

   You are here: :doc:`/dev/index` › :doc:`/dev/one_time_step` ›
   :doc:`/dev/ppc_order` › Grabbing particles from neighbours

**In one sentence:** a particle that leaves its :term:`patch` for a same-level
neighbour is first kept by the patch it left, as a :term:`patch ghost particle`,
so that this patch deposits it; then the neighbour *grabs* it into its domain.

The mechanism
-------------

.. plot::
   :include-source: False
   :caption: Patches A and B of one level, 2 ghost cells (orders 2-3). Follow
      the particle that crosses from A to B: kept and deposited by A, then
      grabbed by B.

   from phare_draw import *
   from matplotlib.patches import Rectangle
   A, B = (0, 0, 6, 4), (6, 0, 12, 4)
   def base(a):
       grid(a, (-2, -2, 13, 6.5), color=C["light"])
       patch(a, A, ghosts=2, label="patch A"); patch(a, B, label="patch B")
       a.add_patch(Rectangle((6, 0), 2, 4, fc=C["E"], alpha=0.3, ec="none", zorder=1))
   fig, (a1, a2, a3) = strip(3, [
       "Push (in A): leavers sorted out of A's domain",
       "A deposits it; moments summed on shared nodes",
       "Grab: B copies it into its domain; A empties the array"], height=2.6)
   base(a1)
   arrow(a1, (4.8, 2.0), (6.7, 2.4), color=C["grey"], lw=1.3, ms=10)
   particles(a1, [(4.8, 2.0)], color=C["ion"], alpha=0.35); particles(a1, [(6.7, 2.4)], color=C["L2"])
   text(a1, 9.6, 2.6, "→ A's patch ghost array", color=C["L2"])
   arrow(a1, (1.0, 3.3), (0.4, 4.7), color=C["grey"], lw=1.3, ms=10)
   text(a1, 0.4, 4.75, "✕", color=C["bad"], bold=True)
   text(a1, 2.2, 5.6, "into the level ghost layer: erased", color=C["bad"], ha="left")
   base(a2)
   particles(a2, [(6.7, 2.4)], color=C["L2"])
   for (x, y) in [(i, j) for i in (6, 7, 8) for j in (1, 2, 3)]:
       mark(a2, "primal", x, y, color=C["rho"], size=7)
   text(a2, 10.4, 2.0, "order 2: 3×3 nodes\non A's grid", color=C["rho"])
   base(a3)
   particles(a3, [(6.7, 2.4)], color=C["L1"])
   text(a3, 9.6, 2.6, "now in B's domain", color=C["L1"])
   text(a3, 3, -1.2, "A's patch ghost array: emptied", color=C["L2"])
   finish(fig, [a1, a2, a3])

* **Push 2 only** (the real push): leavers in a neighbour's cells go to the
  patch ghost array, all leavers leave the domain array
  (``IonUpdater::updateAndDepositAll_``,
  :src:`src/core/numerics/ion_updater/ion_updater.hpp#L258-L289`).
* Domain and patch ghost particles are deposited (L306-L309); density and flux
  are then summed over patch borders (:doc:`/dev/moments`).
* ``fillIonGhostParticles`` runs the grab schedule, then clears every patch
  ghost array (:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L375-L383`).
* Push 1 also exchanges moments (border sums, level ghosts) but **no particles**
  (:src:`src/amr/solvers/solver_ppc.hpp#L620-L664`).
* The receiving patch does not deposit the particle in this step: A already did.

The tricks
----------

Deposit before handing over
^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** A crossing particle must be counted once, though it changes owner.
**Fix.** A deposits it from its patch ghost array; the border sum adds that to
B's copy of the shared nodes (:src:`src/core/numerics/ion_updater/ion_updater.hpp#L306-L309`).

The exchange region
^^^^^^^^^^^^^^^^^^^

**Problem.** B must take every leaver of A that lands in its domain, and only those.
**Fix.** Region = A's box grown by the particle ghost width ∩ B's domain
(``ParticleDomainFromGhostFillPattern``,
:src:`src/amr/data/particles/particles_variable_fill_pattern.hpp#L80-L90`);
copy (same rank) or pack/unpack (other rank) from A's patch ghosts into B's domain
(:src:`src/amr/data/particles/particles_data.hpp#L709-L793`). :pr:`1334`.

One owner per cell
^^^^^^^^^^^^^^^^^^

**Problem.** SAMRAI can make same-level patches overlap: particles would be pushed
and deposited twice.
**Fix.** A cell belongs to the covering patch with the smallest global id
(``makeNonOwnedBoxes``, :src:`src/amr/resources_manager/amr_utils.hpp#L246-L280`);
after the grab each patch drops particles in cells it does not own
(:src:`src/amr/solvers/solver_ppc.hpp#L651-L656`). :pr:`1334`.

Periodic neighbours
^^^^^^^^^^^^^^^^^^^

SAMRAI gives the exchange a shift by the domain size; the particle cell index is
shifted by it (:src:`src/amr/data/particles/particles_data.hpp#L713-L730`).

How far a particle can move
^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** Everything above assumes leavers stay near the patch they left.
**Fix.** Boris stops the run if the new offset ``delta + v Δt/(2Δx)`` exceeds 2
in a half push (:src:`src/core/numerics/pusher/boris.hpp#L157-L172`): at most
about 2 cells forward, 3 backward. The message is only in ``.log/<rank>.out``.
Cell maps cover that margin (``particleMapMargin``,
:src:`src/core/data/grid/gridlayout.hpp#L59-L68`). :pr:`1592`.

.. admonition:: Known issue
   :class: note

   Same-level neighbours are searched with a width of 1 cell
   (:src:`src/amr/resources_manager/amr_utils.hpp#L216-L220`) while the particle
   ghost width is 2 at orders 2-3: a neighbour touching the ghost box only at a
   corner can be missed, and a particle pushed there is erased. Fix proposed in
   :pr:`1593` (corner neighbour search).

Go deeper
---------

* :doc:`/dev/moments` — the border sum that completes the shared nodes.

Related
-------

* :doc:`/dev/level_ghost_particles` · :doc:`/dev/ppc_order` ·
  :doc:`/dev/particles_ions` · :doc:`/dev/ghost_cells`
