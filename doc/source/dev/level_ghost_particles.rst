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
  :src:`src/amr/messengers/refiner.hpp#L161-L166`).
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
  advanced (:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L482-L505`).
* Every fine substep pushes the copy; entrants become domain particles, leavers
  are erased (``IonUpdater::updateAndDepositAll_``,
  :src:`src/core/numerics/ion_updater/ion_updater.hpp#L294-L306`).
* ``fillIonPopMomentGhosts`` deposits Old × (1−α) and New × α of the
  particles the patch *owns*, before the border sums, with
  α = (t_fine − t_start)/(t_end − t_start)
  (:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L418-L455`,
  ``timeInterpCoef_`` :src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L935-L939`).
* ``lastStep``: swap, empty New, reset the copy
  (:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L516-L539`).
* 4 substeps: the fine dt is the coarse dt / 2² (``getMaxFinerLevelDt``,
  :src:`src/amr/multiphysics_integrator.hpp#L420-L428`).

The tricks
----------

Three arrays, not one
^^^^^^^^^^^^^^^^^^^^^

**Problem.** The pushed copy empties during the substeps: nothing refills the
layer from outside, so it cannot give the border moments.
**Fix.** Old and New give the moments (time interpolation); the copy gives
the entrants (:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L442-L452`).

New is split once per coarse step
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** Splitting is a SAMRAI refine :term:`schedule` across levels and
:term:`MPI` ranks: too costly at every substep.
**Fix.** Split once at ``firstStep`` (:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L497`), then interpolate in time; α outside
[0, 1] throws (:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L423-L431`).

One patch deposits each level ghost particle
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. plot::
   :include-source: False
   :caption: The level border is the thick line: refined level below, coarse
      region above. The 2-cell ghost boxes of A and
      B overlap on a band 4 cells wide in the level ghost layer (hatched): both patches hold
      a copy of the particles there. A comes first (its box starts further
      left), so A deposits them; B cedes those cells.

   from phare_draw import *
   from matplotlib.patches import Rectangle
   fig, (a,) = strip(1, ["Two patches below the level border"], height=3, phone=True)
   grid(a, (1, 3, 15, 8), color=C["light"])
   patch(a, (1, 3, 8, 6), color=C["L1"]); patch(a, (8, 3, 15, 6), color=C["L2"])
   text(a, 4.5, 4.5, "A", color=C["L1"], bold=True, fs=14); text(a, 11.5, 4.5, "B", color=C["L2"], bold=True, fs=14)
   a.plot([1, 15], [6, 6], color=C["ink"], lw=2.5, zorder=4)
   a.add_patch(Rectangle((1, 6), 5, 2, fc=C["L1"], alpha=0.25, ec="none", zorder=0.5))
   a.add_patch(Rectangle((10, 6), 5, 2, fc=C["L2"], alpha=0.25, ec="none", zorder=0.5))
   a.add_patch(Rectangle((6, 6), 4, 2, fc="none", ec=C["ink"], hatch="//", lw=1.6, zorder=0.7))
   particles(a, [(7.4, 7.3), (9.2, 6.6)], color=C["ink"])
   text(a, 3.5, 7, "A only", color=C["L1"], bold=True)
   text(a, 12.5, 7, "B only", color=C["L2"], bold=True)
   text(a, 8, 9.0, "both: A deposits", color=C["ink"], bold=True)
   arrow(a, (8, 8.7), (8, 8.1), color=C["ink"], lw=1.4, ms=10)
   finish(fig, [a])

**Problem.** Level ghost particles are split into the whole ghost box of every
patch that covers them, so near a patch border two patches (or more) hold a
copy of the same particle. Each depositing its copy would count it twice once
the border sums add the deposits.
**Fix.** All copies are pushed, but only the *owner* deposits: the first of
the holding patches in lexicographic order of their boxes (lower corner, then
upper corner; ``lexicographicLess``,
:src:`src/core/utilities/box/box.hpp#L159-L167`), identical boxes going to the
smaller patch id. A patch *cedes* the part of its ghost box that lies in the
ghost box of a patch before it (``LevelGhostOwnership``,
``makeLevelGhostOwnership``,
:src:`src/core/numerics/moments/level_ghost_ownership.hpp#L18-L53`);
``depositOwned`` filters on the particle's cell
(:src:`src/core/numerics/moments/moments.hpp#L25-L34`).

* The order is on boxes: a common shift keeps it, so periodic images agree
  on the owner.
* Neighbours searched within 2 × the particle ghost width: two ghost boxes of
  width w overlap only if the patches are less than 2w cells apart
  (``levelGhostNeighbourWidth``,
  :src:`src/core/numerics/moments/level_ghost_ownership.hpp#L32-L38`).
* Built for every refined level with its same-level schedules, so again each
  time a level is replaced (``buildSameLevelSchedules_``,
  :src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L653-L670`);
  debug builds check that every level ghost cell has exactly one owner
  (``checkLevelGhostOwnership``, :src:`src/amr/resources_manager/amr_utils.hpp#L394-L417`).
* The owner deposits *before* the border sums, which add its deposit to the
  other patches (:doc:`/dev/moments`). :pr:`1580` (deposit once, before sums).

Entering a neighbour's domain: erased
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. plot::
   :include-source: False
   :caption: A re-entrant corner of the level border: B is shorter than A.
      The particle's cell is a level ghost cell of A and of B, so both hold a
      copy. Pushed into B's domain, B's copy enters B; A's copy is erased.

   from phare_draw import *
   from matplotlib.patches import Rectangle
   fig, (a1, a2) = strip(2, ["Before the push: a copy on A and on B",
                             "After: B takes it in, A drops its copy"], height=3.0)
   for a in (a1, a2):
       grid(a, (2, 0, 14, 8), color=C["light"])
       a.add_patch(Rectangle((2, 6), 6, 2, fc=C["ghost"], ec="none", zorder=0.4))
       a.add_patch(Rectangle((8, 4), 6, 4, fc=C["ghost"], ec="none", zorder=0.4))
       patch(a, (2, 0, 8, 6), color=C["L1"]); patch(a, (8, 0, 14, 4), color=C["L2"])
       a.plot([2, 8, 8, 14], [6, 6, 4, 4], color=C["ink"], lw=2.5, zorder=4)
       text(a, 5, 3, "A", color=C["L1"], bold=True, fs=14); text(a, 12.5, 1, "B", color=C["L2"], bold=True, fs=14)
       text(a, 11.5, 6.5, "outside the level", color=C["ghost_edge"])
   particles(a1, [(8.6, 4.9)], color=C["L1"], size=40); particles(a1, [(8.9, 4.6)], color=C["L2"], size=40)
   text(a1, 9.3, 5.3, "copies on A and B", ha="left")
   particles(a2, [(8.9, 3.3)], color=C["L2"], size=40)
   text(a2, 9.3, 3.3, "B's copy: enters, kept", ha="left", color=C["L2"], bold=True)
   text(a2, 8.45, 3.3, "✕", color=C["bad"], bold=True, fs=14, halo=False)
   text(a2, 8.3, 2.4, "A's copy: erased", ha="left", color=C["bad"], bold=True)
   arrow(a2, (8.9, 4.8), (8.9, 3.6), color=C["grey"], lw=1.2, ms=9)
   finish(fig, [a1, a2])

**Problem.** After the push, A kept the level ghost particles still in its
ghost box and outside its own domain. At a re-entrant corner of the level
border, such a particle can step into the domain of the neighbour B, which
holds its own copy and takes it in. A kept its copy in its level ghost array:
when the particle later moved into A's domain, A took that copy in while B's
copy also reached A through the patch exchange, one particle twice.
**Fix.** A keeps only the copies in its level ghost layer: its ghost box minus
its domain *and* minus the neighbours' domains (``inLevelGhostLayer``,
``nonLevelGhostBox``,
:src:`src/core/numerics/ion_updater/ion_updater.hpp#L167-L174`;
``makeNonLevelGhostBoxFor``,
:src:`src/amr/resources_manager/amr_utils.hpp#L233-L251`). A's copy is neither
taken into A's domain nor kept: it is erased
(:src:`src/core/numerics/ion_updater/ion_updater.hpp#L294-L306`).
:pr:`1579` (corner duplicate fix).

Domain particles that leave the level are dropped
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** Keeping them would count that plasma twice (also in Old/New).
**Fix.** Only leavers in the *non-level-ghost box* become patch ghosts; the
others are erased (:src:`src/core/numerics/ion_updater/ion_updater.hpp#L268-L292`).

Entrants of the first push are deposited, not kept
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** The first push only predicts moments; particles stay at time n
(:doc:`/dev/ppc_order`).
**Fix.** A temporary copy is pushed; entrants are deposited and thrown away
(:src:`src/core/numerics/ion_updater/ion_updater.hpp#L213-L236`).

Regrid fills Old by splitting only
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** A SAMRAI regrid schedule for Old caused occasional MPI failures.
**Fix.** ``regrid`` uses a plain fill (``fillTransient``)
(:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L268-L282`).
:issue:`604` (levelGhostOld regrid failures).

Restarts keep them
^^^^^^^^^^^^^^^^^^

Old, New and the copy are written to restart files; patch ghosts are not
(:src:`src/amr/data/particles/particles_data.hpp#L137-L140`).

A limit: fast particles
^^^^^^^^^^^^^^^^^^^^^^^

The copy is filled once per coarse step
(:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L536`), so
a particle farther out than the layer cannot enter during that coarse step.
Entry is complete only if particles move less than the layer width over 4 fine
steps. Border moments are not affected (they come from Old and New).

Go deeper
---------

* :doc:`/dev/particle_splitting` — how coarse particles become fine ones.

Related
-------

* :doc:`/dev/one_time_step` · :doc:`/dev/ppc_order` ·
  :doc:`/dev/grabbing_particles` · :doc:`/dev/moments` · :doc:`/dev/ghost_cells`
