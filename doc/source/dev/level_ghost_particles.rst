.. _dev-level-ghost-particles:

=====================
Level ghost particles
=====================

.. container:: breadcrumb

   You are here: :doc:`/dev/index` › :doc:`/dev/one_time_step` › Level ghost particles

**In one sentence:** level ghost particles are coarse particles, split onto the
fine grid and kept in a thin layer just outside a refined :term:`level`, so that
plasma can flow *into* the level and the ion :term:`moments` near its border are
complete.

Where they live
---------------

.. plot::
   :include-source: False
   :caption: Patch A's particle ghost box (orange, dashed outline) is cut in two: the part over
      patch B is filled by B (patch ghosts, blue); the part outside the level is
      the level ghost layer, filled from the coarser level.

   from phare_draw import *
   from matplotlib.patches import Rectangle
   fig, (a,) = strip(1, ["A refined level made of patches A and B"], height=3.4)
   grid(a, (4, 4, 28, 20), color=C["light"])
   patch(a, (8, 8, 16, 16), ghosts=2, label="patch A")
   patch(a, (16, 8, 24, 16), label="patch B")
   a.add_patch(Rectangle((16, 8), 2, 8, fc=C["E"], alpha=0.35, ec="none", zorder=1))
   text(a, 6.5, 21.3, "level ghost layer of A", color=C["ghost_edge"], bold=True)
   text(a, 22, 5.0, "patch ghosts of A (from B)", color=C["E"], bold=True)
   arrow(a, (6.5, 20.7), (7.0, 16.5), color=C["ghost_edge"], lw=1.5)
   arrow(a, (20, 5.6), (17, 9), color=C["E"], lw=1.5)
   finish(fig, [a])

* Width: the particle ghost width, 1 fine cell at interpolation order 1, 2 at
  orders 2 and 3 (``particle_ghost_width``,
  :src:`src/core/models/options/hybrid_options.hpp#L25`).
* Filled with SAMRAI's ``PatchLevelBorderFillPattern``: patch box grown by the
  ghost width, minus the whole level (``RefinerType::LevelBorderParticles``,
  :src:`src/amr/messengers/refiner.hpp#L188-L193`).
* Root level (level 0) has none: its border is the physical boundary.

Why they are needed
-------------------

.. plot::
   :include-source: False
   :caption: 1D, order 1. The node on the level border (red) gets a share of
      density from a particle inside the level (blue) and one outside it
      (orange). Without the outside particle, the border density is too low.

   import numpy as np
   from phare_draw import *
   fig, (a,) = strip(1, ["Two particle clouds reaching the border node"], height=2.6, aspect=False)
   x = np.linspace(0, 9, 901)
   hat = lambda c: np.clip(1 - np.abs(x - c), 0, None) * 1.6
   a.axvspan(5, 9, color=C["ghost"], lw=0)
   a.axvline(5, color=C["ink"], lw=2)
   a.fill_between(x, hat(4.4), color=C["E"], alpha=0.35)
   a.fill_between(x, hat(5.5), color=C["ghost_edge"], alpha=0.6)
   for n in range(10):
       mark(a, "primal", n, 0, color=C["bad"] if n == 5 else C["grey"], size=11 if n == 5 else 7)
   text(a, 2.5, 1.9, "refined level"); text(a, 7, 1.9, "outside (coarse)")
   finish(fig, [a])

* No particle ever enters the level from the coarse region around it.
* Nodes near the border miss the share of the clouds sitting outside.

One coarse step, three arrays
-----------------------------

.. plot::
   :include-source: False
   :caption: The three level ghost arrays of one patch over one coarse step
      (4 fine substeps). Old and New are never pushed; the copy is.

   from phare_draw import *
   from matplotlib.patches import Rectangle
   rows = {"Old": 2.5, "New": 1.5, "copy": 0.5}
   old = [8.3, 8.9, 9.4, 9.8]; new = [8.6, 9.1, 9.5]
   def base(a):
       a.add_patch(Rectangle((8, 0), 2, 3, fc=C["ghost"], ec="none", zorder=0))
       grid(a, (4, 0, 12, 3), color=C["light"])
       a.plot([8, 8], [0, 3], color=C["ink"], lw=2.5, zorder=3)
       for name, y in rows.items():
           text(a, 4.1, y, name, ha="left", bold=True)
   fig, (a1, a2, a3) = strip(3, [
       "First substep: Old split at t, New at t+Δt",
       "Substep k: push the copy; moments (1−α)·Old + α·New",
       "Last substep: Old ← New, New emptied, copy ← Old"], height=2.0)
   base(a1)
   particles(a1, [(v, 2.5) for v in old], color=C["grey"])
   particles(a1, [(v, 1.5) for v in new], color=C["L1"])
   particles(a1, [(v, 0.5) for v in old], color=C["L2"])
   text(a1, 6, 3.4, "fine level"); text(a1, 9, 3.4, "level ghost layer", color=C["ghost_edge"])
   base(a2)
   particles(a2, [(v, 2.5) for v in old], color=C["grey"])
   particles(a2, [(v, 1.5) for v in new], color=C["L1"])
   text(a2, 10.6, 2.5, "×(1−α)", ha="left"); text(a2, 10.6, 1.5, "×α", ha="left")
   moved = [(8.3, 7.4), (8.9, 8.5), (9.4, 9.0), (9.8, 10.4)]
   for p, q in moved:
       arrow(a2, (p, 0.5), (q, 0.5), color=C["grey"], lw=1.2, ms=9)
   particles(a2, [(7.4, 0.5)], color=C["ion"]); text(a2, 6.3, 0.5, "enters: domain", color=C["ion"])
   particles(a2, [(8.5, 0.5), (9.0, 0.5)], color=C["L2"])
   text(a2, 10.4, 0.5, "✕", color=C["bad"], bold=True); text(a2, 10.6, 0.15, "left: erased", ha="left", color=C["bad"])
   base(a3)
   particles(a3, [(v, 2.5) for v in new], color=C["grey"])
   text(a3, 9, 1.5, "(empty)", color=C["grey"])
   particles(a3, [(v, 0.5) for v in new], color=C["L2"])
   finish(fig, [a1, a2, a3])

* ``firstStep``: New is split from the coarse level, which has already
  advanced (:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L508-L532`).
* Every fine substep pushes the copy; entrants become domain particles, leavers
  are erased (``IonUpdater::updateAndDepositAll_``,
  :src:`src/core/numerics/ion_updater/ion_updater.hpp#L291-L304`).
* ``fillIonPopMomentGhosts`` deposits Old × (1−α) and New × α, with
  α = (t_fine − t_start)/(t_end − t_start)
  (:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L457-L495`,
  ``timeInterpCoef_`` :src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L949-L953`).
* ``lastStep``: swap, empty New, reset the copy
  (:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L542-L565`).
* 4 substeps: the fine dt is the coarse dt / 2² (``getMaxFinerLevelDt``,
  :src:`src/amr/multiphysics_integrator.hpp#L420-L428`).

The tricks
----------

Three arrays, not one
^^^^^^^^^^^^^^^^^^^^^

**Problem.** The pushed copy empties during the substeps: nothing refills the
layer from outside, so it cannot give the border moments.
**Fix.** Old and New give the moments (time interpolation); the copy gives
the entrants (:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L481-L492`).

New is split once per coarse step
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** Splitting is a SAMRAI refine :term:`schedule` across levels and
:term:`MPI` ranks: too costly at every substep.
**Fix.** Split once at ``firstStep`` (:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L523`), then interpolate in time; α outside
[0, 1] throws (:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L462-L470`).

Domain particles that leave the level are dropped
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** Keeping them would count that plasma twice (also in Old/New).
**Fix.** Only leavers in the *non-level-ghost box* become patch ghosts; the
others are erased (:src:`src/core/numerics/ion_updater/ion_updater.hpp#L265-L289`).

Entrants of the first push are deposited, not kept
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** The first push only predicts moments; particles stay at time n
(:doc:`/dev/ppc_order`).
**Fix.** A temporary copy is pushed; entrants are deposited and thrown away
(:src:`src/core/numerics/ion_updater/ion_updater.hpp#L208-L233`).

Regrid fills Old by splitting only
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** A SAMRAI regrid schedule for Old caused occasional MPI failures.
**Fix.** ``regrid`` uses a plain fill (``fillTransient``)
(:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L266-L280`).
:issue:`604` (levelGhostOld regrid failures).

Restarts keep them
^^^^^^^^^^^^^^^^^^

Old, New and the copy are written to restart files; patch ghosts are not
(:src:`src/amr/data/particles/particles_data.hpp#L137-L140`).

A limit: fast particles
^^^^^^^^^^^^^^^^^^^^^^^

The copy is filled once per coarse step
(:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L560-L562`), so
a particle farther out than the layer cannot enter during that coarse step.
Entry is complete only if particles move less than the layer width over 4 fine
steps. Border moments are not affected (they come from Old and New).

.. admonition:: Known issue
   :class: note

   Near a corner of the level border, two neighbour patches both hold a copy of
   the same level ghost particle; it can enter twice. Fix proposed in
   :pr:`1579` (corner duplicate fix).

Go deeper
---------

* :doc:`/dev/particle_splitting` — how coarse particles become fine ones.

Related
-------

* :doc:`/dev/one_time_step` · :doc:`/dev/ppc_order` ·
  :doc:`/dev/grabbing_particles` · :doc:`/dev/moments` · :doc:`/dev/ghost_cells`
