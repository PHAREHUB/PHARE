.. _dev-level-ghost-particles:

=====================
Level ghost particles
=====================

.. container:: breadcrumb

   You are here: :doc:`/dev/index` › :doc:`/dev/one_time_step` › Level ghost particles

**In one sentence:** level ghost particles are copies of coarse-level particles,
split onto the fine grid and placed in a thin layer just outside a refined level,
so that particles can flow *into* the refined level and the ion moments near its
border are complete.

Why they exist
--------------

A refined level covers only part of the domain. Its border is not a physical
boundary: plasma keeps flowing across it. Two things go wrong at that border if
nothing is done:

* **Nothing comes in.** Fine particles that leave the level are simply lost, and
  no particle ever enters from the coarse region around it.
* **Moments are wrong near the border.** The *ion moments* are the density and
  the flux of the ions on the grid nodes. They are computed by *depositing* the
  particles: each particle is a small cloud (its *shape*, one to a few cells wide
  depending on the interpolation order chosen by the user), and every node adds
  up the share of each cloud that covers it (:doc:`/dev/moments`). Nodes near the
  border are also covered by clouds sitting *outside* the level. Without those
  particles, the density drops at the border:

  .. plot::
     :include-source: False
     :caption: 1D, interpolation order 1 (each cloud is a triangle two cells
        wide). The node on the level border (red) gets a share from a
        particle inside the level (blue) and one from a particle outside it
        (orange). Without the outside particle, that node misses part of its
        density. With wider shapes (orders 2 and 3), outside particles reach
        further into the level.

     import numpy as np
     import matplotlib.pyplot as plt

     fig, ax = plt.subplots(figsize=(6, 3))
     x = np.linspace(0, 9, 901)
     hat = lambda c: np.clip(1 - np.abs(x - c), 0, None)
     ax.axvspan(5, 9, color="#f4a259", alpha=0.25, lw=0)
     ax.axvline(5, color="k", lw=2)
     ax.fill_between(x, hat(4.4), color="#2a7ab9", alpha=0.35)
     ax.fill_between(x, hat(5.5), color="#f4a259", alpha=0.8)
     for n in range(0, 10):
         ax.plot([n], [0], "o", color="red" if n == 5 else "0.4",
                 ms=10 if n == 5 else 6, zorder=5, clip_on=False)
     ax.text(2.5, 1.12, "refined level", ha="center", fontsize=12)
     ax.text(7, 1.12, "outside (coarse)", ha="center", fontsize=12)
     ax.text(5.1, -0.32, "level border", fontsize=11)
     ax.set_ylim(-0.4, 1.3); ax.set_xlim(0, 9)
     ax.set_xticks(range(10)); ax.set_yticks([])
     ax.set_xlabel("grid nodes", fontsize=12)
     ax.tick_params(labelsize=11)
     for side in ("top", "right", "left"):
         ax.spines[side].set_visible(False)
     fig.tight_layout()

The coarser level does have particles there. PHARE takes them, splits each one
into several finer particles (:doc:`/dev/particle_splitting`), and keeps the
pieces that fall in a layer one or two fine cells wide around the level. Those
are the level ghost particles.

The picture
-----------

Where they live: a refined level made of two patches, A and B, seen from patch A.

.. plot::
   :include-source: False
   :caption: Patch A's particle ghost box (dashed) is cut in two. The part over
      patch B is a *patch ghost* region, filled by the same-level neighbour.
      The part outside the level (orange) is the *level ghost* layer, filled from
      the coarser level. Width: 2 fine cells here (interpolation order 2 or 3;
      1 cell at order 1).

   import matplotlib.pyplot as plt
   from matplotlib.patches import Rectangle

   fig, ax = plt.subplots(figsize=(6, 5))
   g = 2                      # particle ghost width (fine cells)
   A = (8, 8, 8, 8)           # x0, y0, width, height in fine cells
   B = (16, 8, 8, 8)
   # coarse grid (2 fine cells per coarse cell)
   for x in range(2, 31, 2):
       ax.axvline(x, color="0.85", lw=0.8, zorder=0)
   for y in range(2, 23, 2):
       ax.axhline(y, color="0.85", lw=0.8, zorder=0)
   # level ghost layer of A = ghost box minus the level (A and B)
   gx0, gy0, gw, gh = A[0] - g, A[1] - g, A[2] + 2 * g, A[3] + 2 * g
   ax.add_patch(Rectangle((gx0, gy0), gw, gh, color="#f4a259", alpha=0.85, zorder=1))
   for P in (A, B):
       ax.add_patch(Rectangle(P[:2], P[2], P[3], color="white", zorder=2))
   # patch ghost part of A's ghost box (over B)
   ax.add_patch(Rectangle((B[0], B[1]), g, B[3], color="#7fb7e6", alpha=0.9, zorder=3))
   # fine grid inside the level
   for P, name in ((A, "patch A"), (B, "patch B")):
       x0, y0, w, h = P
       for x in range(x0, x0 + w + 1):
           ax.plot([x, x], [y0, y0 + h], color="0.75", lw=0.5, zorder=4)
       for y in range(y0, y0 + h + 1):
           ax.plot([x0, x0 + w], [y, y], color="0.75", lw=0.5, zorder=4)
       ax.add_patch(Rectangle(P[:2], w, h, fill=False, lw=2, color="k", zorder=5))
       ax.text(x0 + w / 2, y0 + h / 2, name, ha="center", va="center", fontsize=13, zorder=6)
   ax.add_patch(Rectangle((gx0, gy0), gw, gh, fill=False, ls="--", lw=1.5,
                          color="#333333", zorder=6))
   ax.annotate("level ghost layer of A\n(from the coarser level)", xy=(7, 12),
               xytext=(0, 23), fontsize=12,
               arrowprops=dict(arrowstyle="->", color="#333333"))
   ax.annotate("patch ghost region of A\n(from patch B)", xy=(17, 12),
               xytext=(15, 1), fontsize=12,
               arrowprops=dict(arrowstyle="->", color="#333333"))
   ax.text(24.5, 18.5, "coarse level\n(grey grid)", fontsize=11, color="0.4")
   ax.set_xlim(0, 30); ax.set_ylim(0, 26)
   ax.set_aspect("equal"); ax.axis("off")
   fig.tight_layout()

What happens to them during one coarse step: three particle arrays, each with
its own job (click a box to go deeper).

.. mermaid::

   flowchart TD
     C["coarse level<br/>domain particles"] -->|split| O["levelGhostParticlesOld<br/>(coarse time t)"]
     C -->|"split, at the 1st fine substep"| N["levelGhostParticlesNew<br/>(coarse time t + Δt)"]
     O -->|copy| P["levelGhostParticles<br/>(pushable copy)"]
     P -->|"pushed each substep"| D["fine domain particles<br/>(those that enter)"]
     O -->|"deposit × (1 − α)"| M["ion moments<br/>near the level border"]
     N -->|"deposit × α"| M
     N -.->|"last substep:<br/>New becomes Old"| O
     click C href "particle_splitting.html" "How coarse particles are split"
     click O href "particle_splitting.html" "How coarse particles are split"
     click N href "particle_splitting.html" "How coarse particles are split"
     click M href "moments.html" "Moments: deposition and sum borders"
     click D href "grabbing_particles.html" "Particles moving between patches"

How it works
------------

A fine level takes several smaller steps for each coarse step. With PHARE's
refinement ratio of 2, the fine time step is the coarse one divided by
:math:`2^2 = 4`, because whistler waves need :math:`\Delta t \propto \Delta x^2`
(``MultiPhysicsIntegrator::getMaxFinerLevelDt``,
:src:`src/amr/multiphysics_integrator.hpp#L420-L428`). So there are 4 fine
substeps per coarse step. The coarse level always advances first, so during the
fine substeps the coarse particles already sit at the *end* of the coarse step.

Each patch of a refined level stores, per ion population, five particle arrays
(``ParticlesData``, :src:`src/amr/data/particles/particles_data.hpp#L460-L466`):
domain, patch ghost, and three for level ghosts:

.. list-table::
   :header-rows: 1
   :widths: 25 20 55

   * - Array
     - Time
     - Job
   * - ``levelGhostParticlesOld``
     - start of the coarse step
     - deposit moments near the border; never pushed
   * - ``levelGhostParticlesNew``
     - end of the coarse step
     - deposit moments near the border; never pushed
   * - ``levelGhostParticles``
     - follows the fine substeps
     - a copy of Old that *is* pushed; the particles that enter the level
       become domain particles

The cycle over one coarse step, for a level above the root:

#. **First fine substep** — ``HybridHybridMessengerStrategy::firstStep``
   (:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L508-L532`)
   splits the coarse domain particles, now at the end of the coarse step, into
   ``levelGhostParticlesNew`` (L523). It also stores the start and end times of
   the coarse step (L528-L529). It is called by
   ``MultiPhysicsIntegrator::advanceLevel`` only on the first substep
   (:src:`src/amr/multiphysics_integrator.hpp#L524-L531`).
#. **Every substep, push** — the hybrid solver pushes particles twice per
   substep: a first push only predicts the moments, the second really moves the
   particles (:doc:`/dev/ppc_order`). In the second push, the pushable copy
   ``levelGhostParticles`` is moved with the fine fields. Particles that land in a cell
   owned by the patch are moved to the domain array; particles that leave the
   ghost layer are erased (``IonUpdater::updateAndDepositAll_``,
   :src:`src/core/numerics/ion_updater/ion_updater.hpp#L291-L304`). This is how
   plasma enters the refined level.
#. **Every substep, moments** — after each push the solver calls
   ``fillIonPopMomentGhosts``
   (:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L457-L495`).
   It deposits Old with weight :math:`1-\alpha` and New with weight
   :math:`\alpha` (L485-L491), where

   .. math::

      \alpha = \frac{t_\text{fine} - t_\text{coarse}^\text{start}}
                    {t_\text{coarse}^\text{end} - t_\text{coarse}^\text{start}}

   is the position of the fine time inside the coarse step
   (``timeInterpCoef_``, L949-L953). This is a linear interpolation in time
   of the coarse plasma at the border.
#. **Last fine substep** — ``lastStep`` (L542-L565) swaps New and Old, empties
   New, and resets the pushable copy to the new Old (L560-L562). The next coarse
   step starts from fresh coarse particles.

With 4 substeps, :math:`\alpha` takes the values 1/4, 1/2, 3/4 and 1.

.. mermaid::

   sequenceDiagram
     participant C as coarse level
     participant F as fine level (one patch)
     C->>C: advance t → t+Δt
     C->>F: firstStep: split coarse particles → New (at t+Δt)
     loop 4 fine substeps
       F->>F: push the pushable copy, entrants join the domain
       F->>F: moments at border = (1−α)·Old + α·New
     end
     F->>F: lastStep: Old ← New, New emptied, pushable ← Old

When a level is created or rebuilt, ``initLevel`` and ``regrid`` fill Old by
splitting and copy it into the pushable array
(:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L336-L341`,
L276-L280). New is only filled at the next first substep.

The tricks
----------

Why three arrays and not one
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** The pushed level ghost particles cannot be used for the moments.
Nothing refills the ghost layer from outside during the substeps: particles that
leave it are erased, and no particle enters it from the coarse side. After a
substep the layer is partly empty, and the density near the border would drop.

**Fix.** Two jobs, two kinds of arrays. Old and New are never pushed: they give
the *moments*, by time interpolation between two complete snapshots of the
coarse plasma. The pushable copy gives the *particles* that enter the level.
The interpolation is done in ``fillIonPopMomentGhosts``
(:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L481-L492`).

New is filled at the first substep, not at every substep
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** Splitting means a parallel communication across the coarse and
fine levels: a SAMRAI *refine schedule*, the plan of which data goes from which
coarse patch to which fine patch, possibly on another MPI rank (one of the
parallel processes of the run). Doing it at every fine substep
would be costly.

**Fix.** The coarse level is ahead in time, so its particles at the end of the
coarse step are already known at the first substep. They are split once into
New (L523); the next substeps only interpolate in time, with no communication.
``fillIonPopMomentGhosts`` throws an error if :math:`\alpha` falls outside
[0, 1] on a refined level (L462-L470): that would mean the times stored by
``firstStep`` do not frame the fine time.

Only the layer outside the level is filled
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** A patch's particle ghost box also covers its same-level
neighbours. Particles there must come from those neighbours (they are fine
particles already), not from the coarse level.

**Fix.** The level ghost refiners use SAMRAI's ``PatchLevelBorderFillPattern``
(``Refiner`` for ``RefinerType::LevelBorderParticles``,
:src:`src/amr/messengers/refiner.hpp#L188-L193`). It grows each patch box by the
ghost width and removes the whole level from it: only the orange region of the
picture is filled. The width is the particle ghost width: 1, 2, 2 fine cells for
interpolation orders 1, 2, 3 (``particle_ghost_width``,
:src:`src/core/models/options/hybrid_options.hpp#L25`).

Domain particles that cross into the level ghost layer are dropped
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** A fine domain particle that moves out of the level is now in the
region the coarse level represents. Keeping it, or depositing it there, would
count that plasma twice: once through it, once through Old and New.

**Fix.** After the push, only domain particles that land in the patch's
*non-level-ghost box* (its domain plus the parts of the ghost box covered by
same-level neighbours) are kept as patch ghosts and deposited; those in the
level ghost layer are erased (``IonUpdater::updateAndDepositAll_``,
:src:`src/core/numerics/ion_updater/ion_updater.hpp#L265-L289`). The same
selection is used in the first, moments-only push (L198-L205).

Entering level ghosts are deposited but not kept in the first push
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** As said above, the first of the two pushes only predicts the
moments; particles must stay at the start time (see :doc:`/dev/ppc_order`).

**Fix.** In the first push (``IonUpdater::updateAndDepositDomain_``,
:src:`src/core/numerics/ion_updater/ion_updater.hpp#L208-L233`), a temporary
copy of the pushable array is pushed, and only the particles that enter owned
domain cells are deposited. They are not copied anywhere: the comment at
L223-L230 explains that the border nodes get their share from Old and New.

Regridding fills Old by splitting only
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** When a level is rebuilt, a SAMRAI *regrid* schedule copies data
from the old level where it overlaps the new one, and refines elsewhere. Used
for Old, it led to occasional failures in SAMRAI's parallel communication.

**Fix.** ``regrid`` fills Old with a plain fill schedule
(``lvlGhostPartOldRefiners_.fillTransient``), so Old always comes from splitting
coarse particles, as for a brand new level
(:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L266-L280`).

**History.** #604 (levelGhostOld regrid failures).

Restart files keep the level ghost arrays
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The three level ghost arrays are written to restart files together with the
domain particles; patch ghost particles are not
(``ParticlesData::putToRestart``,
:src:`src/amr/data/particles/particles_data.hpp#L137-L140`). A restarted run
therefore resumes with the same Old and pushable particles.

.. admonition:: Known issue
   :class: note

   Near a corner of the level border, the ghost boxes of two neighbouring
   patches P and Q overlap outside the level, so both hold a copy of the same
   level ghost particle. When it moves into Q's domain, Q's copy becomes a
   domain particle, but P's copy is still in P's ghost box and is kept. Later
   it can reach P's domain and be handed to Q again: the particle is counted
   twice. It was reported in 2D at interpolation order 1. A fix is proposed in
   `#1579 <https://github.com/PHAREHUB/PHARE/pull/1579>`_ (corner duplicate fix).

Go deeper
---------

* :doc:`/dev/particle_splitting` — how one coarse particle becomes several fine
  ones, and how only those in the layer are kept.

Related
-------

* :doc:`/dev/one_time_step` — where ``firstStep``, the pushes and ``lastStep``
  sit in a time step.
* :doc:`/dev/ppc_order` — the order of pushes and depositions in the solver.
* :doc:`/dev/grabbing_particles` — the same-level counterpart: patch ghost
  particles.
* :doc:`/dev/moments` — deposition and the sums over patch borders.
* :doc:`/dev/ghost_cells` — the same split between patch and level ghosts, for
  fields.
