.. _dev-messengers:

=========================
Messengers and schedules
=========================

.. container:: breadcrumb

   You are here: :doc:`/dev/index` › :doc:`/dev/one_time_step` › Messengers

**In one sentence:** a :term:`messenger` moves data between the patches of a
level and between a level and the coarser one; it hides SAMRAI's algorithms,
:term:`schedules <schedule>` and :term:`fill patterns <fill pattern>` behind
calls such as ``fillMagneticGhosts(B, level, time)``.


Which messenger
---------------

.. list-table::
   :header-rows: 1
   :widths: 22 22 56

   * - Coarser level
     - Level
     - Messenger
   * - hybrid
     - hybrid
     - ``HybridMessenger`` with ``HybridHybridMessengerStrategy``
   * - MHD
     - MHD
     - ``MHDMessenger``
   * - MHD
     - hybrid
     - ``HybridMessenger`` with ``MHDHybridMessengerStrategy``: all methods
       empty at this commit (:src:`src/amr/messengers/mhd_hybrid_messenger_strategy.hpp`)

* One object per model pair, shared by all levels using it; level 0 uses
  (its model, its model) (``registerMessengers_``,
  :src:`src/amr/multiphysics_integrator.hpp#L820-L861`).
* Created once at setup by the ``Simulator`` (``registerAndSetupMessengers``,
  :src:`src/simulator/simulator.hpp#L331`).


What a schedule is
------------------

.. plot::
   :include-source: false
   :alt: Two level-1 patches over level 0; arrows show the transfers a ghost-fill schedule does for patch A.

   from phare_draw import *

   fig, (ax,) = strip(1, ["the transfers of one ghost fill, for patch A"])
   grid(ax, (0, 0, 16, 7), color=C["light"])
   patch(ax, (4, 2, 8, 5), ghosts=1, label="patch A")
   patch(ax, (8, 2, 12, 5), label="patch B")
   text(ax, 14.2, 6.4, "level 0", color=C["L0"])
   arrow(ax, (9.6, 3.5), (8.5, 3.5), color=C["E"], lw=2.2)
   text(ax, 12.6, 1.0, "copy from patch B (same level,\nmaybe another rank)", color=C["E"])
   arrow(ax, (11.0, 1.6), (9.2, 3.2), color=C["E"], lw=1.2)
   arrow(ax, (2.0, 6.4), (3.6, 4.6), color=C["L2"], lw=2.2)
   text(ax, 1.6, 6.75, "refine from level 0", color=C["L2"])
   ax.plot([-0.3, 16.6], [0.2, 7.4], alpha=0)
   finish(fig, [ax])

Look at patch A's ghost layer: where patch B is, values are copied; elsewhere
they are refined from level 0. A schedule is this list of transfers, for these
patches, on every MPI :term:`rank`.

* **Algorithm** (what to move): quantities, refine operator, fill pattern;
  built once at setup (``registerQuantities``,
  :src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L166-L210`).
* **Schedule** (who sends what to whom): an algorithm applied to a level;
  costly and collective to build; valid while its levels exist.
* **Fill**: ``schedule->fillData(time)``, run by the messenger at each call.
* PHARE wraps them in ``Refiner`` / ``RefinerPool`` (and ``Synchronizer`` for
  coarsening) (:src:`src/amr/messengers/refiner.hpp`,
  :src:`src/amr/messengers/synchronizer.hpp`).


The kinds of transfer
---------------------

``RefinerType`` (:src:`src/amr/messengers/refiner.hpp#L20-L32`) decides how a
schedule is built (``makeSchedule_``, :src:`src/amr/messengers/refiner.hpp#L71-L178`):

.. list-table::
   :header-rows: 1
   :widths: 30 40 30

   * - Type
     - Schedule
     - Hybrid refiners
   * - ``GhostField``
     - copy from the level, refine elsewhere (:doc:`/dev/ghost_cells`)
     - ``magGhostsRefiners_``, ``elecGhostsRefiners_``, ``momentGhostsRefiners_``
       (ion charge density and E-edge current, refined levels, :doc:`/dev/moments`)
   * - ``PatchGhostField``
     - copy from the level only
     - ``chargeDensityPatchGhostsRefiners_``, ``chargeFluxEPatchGhostsRefiners_``
   * - ``InitField``
     - refine everything from the coarser level
     - ``electricInitRefiners_``
   * - ``InitInteriorPart``
     - split coarse particles into the new level's cells
     - ``domainParticlesRefiners_``
   * - ``LevelBorderParticles``
     - split coarse particles just outside the level
       (:doc:`/dev/level_ghost_particles`)
     - ``lvlGhostPartOldRefiners_``, ``lvlGhostPartNewRefiners_``
   * - ``ExteriorGhostParticles``
     - same-level particle exchange (:doc:`/dev/grabbing_particles`)
     - ``domainGhostPartRefiners_``
   * - ``PatchFieldBorderSum`` (``Vec``, ``Tensor``)
     - add neighbours' values on shared ghost-box nodes, in place, in one order
       (:doc:`/dev/moments`)
     - ``popDensityBorderSumRefiners_``, ``popFluxEBorderSumRefiners_``; flux and
       momentum tensor for outputs in ``HybridLazyMoments`` (:doc:`/dev/output_moments`)
   * - ``PatchVecFieldBorderMax``
     - same nodes, keep the larger value
     - ``magMaxRootRefiners_`` (level-0 B, once at init)
   * - ``LevelBorderField``
     - level ghosts only; unused
     - (none)

The hybrid messenger's calls during a step (``prepareStep``, ``firstStep``,
``fill...``, ``lastStep``, ``synchronize``, ``reflux``, ``postSynchronize``)
are listed in order on :doc:`/dev/one_time_step` and
:ref:`dev-one-time-step-sync`.


.. _dev-messengers-schedules:

When schedules are built
------------------------

.. list-table::
   :header-rows: 1
   :widths: 22 38 40

   * - Kind
     - Examples
     - Built
   * - same-level
     - patch ghosts, moment sums, particle exchange
     - when level L is created; rebuilt if L was replaced
       (:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L699-L743`)
   * - inter-level
     - B, E and ion moment ghosts, level ghost particles, coarsening, reflux
     - in ``onHierarchyChanged``; rebuilt if any level 0..L was replaced
       (:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L747-L765`)
   * - one fill only
     - filling a new level; sums of the output moments
     - at fill time, then dropped (``fillTransient``,
       :src:`src/amr/messengers/refiner.hpp#L192-L202`)

Calls: ``registerLevel`` / ``onHierarchyChanged``
(:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L215-L250`),
from ``initializeLevelData`` (:src:`src/amr/multiphysics_integrator.hpp#L328-L331`)
and ``resetHierarchyConfiguration`` (:src:`src/amr/multiphysics_integrator.hpp#L384-L420`).


The tricks
----------

Rebuilt exactly when their levels change
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** After a level-0 rebalance SAMRAI only reports level 0, but level
2 schedules also read level 0; rebuilding everything every regrid is costly.

**Fix.** Each schedule set records the levels it was built on (weak
pointers); ``onHierarchyChanged`` runs for every level and rebuilds only stale
sets (``builtOnCurrentLevels_``,
:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L770-L789`;
:src:`src/amr/multiphysics_integrator.hpp#L414-L420`). :pr:`1419`

Level 0 ghost refiners are same-level
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** B and E ghost refiners are inter-level, but level 0 has no
coarser level.

**Fix.** On the root level they are registered as same-level schedules: copy
plus boundary conditions
(:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L726-L732`).

Init schedules read only the coarser level
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** A new level has no data to copy from its neighbours.

**Fix.** ``InitField`` uses no same-level source (``nullptr``)
(:src:`src/amr/messengers/refiner.hpp#L117-L126`); new particles go in cells
only, never in ghosts, so ghosts stay exact copies of neighbours
(:src:`src/amr/messengers/refiner.hpp#L130-L143`).

Ghosts after init and regrid
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** Init/regrid schedules read coarser levels at their latest time,
which can be ahead of a level created mid-step.

**Fix.** MHD: ``fillGhostsAtInitTime_`` / ``refillGhostsAtInitTime_``
refill them with time interpolation
(:src:`src/amr/messengers/mhd_messenger.hpp#L462-L489`). :pr:`1494`. Hybrid: only the
ion moment ghosts Ohm reads are refilled, with a one-fill schedule, the ions stamped
with the init time so the time interpolation uses it (``fillIonMomentGhosts``,
:src:`src/amr/level_initializer/hybrid_level_initializer.hpp#L123-L132`); B and E ghosts
are next filled in the first sub-step.


Go deeper
---------

* :doc:`/dev/ghost_cells`: what the ghost fills fill.
* :doc:`/dev/fill_patterns`: which nodes each transfer writes.
* :doc:`/dev/field_refinement`: the operators that compute fine values.
* :doc:`/dev/coarsening`: the synchronizers.
* :doc:`/dev/level_ghost_particles` and :doc:`/dev/grabbing_particles`: the particle transfers.

Related
-------

* :doc:`/dev/one_time_step`: when each messenger call happens.
* :doc:`/dev/tagging_regridding`: what creates and replaces levels.
* :doc:`/dev/load_balancing`: level-0 rebalancing.
