.. _dev-messengers:

=========================
Messengers and schedules
=========================

.. container:: breadcrumb

   You are here: :doc:`/dev/index` › :doc:`/dev/one_time_step` › Messengers

**In one sentence:** a :term:`messenger` is the object that moves data between
the patches of a level and between a level and the next coarser one; it hides
SAMRAI's machinery (algorithms, :term:`schedules <schedule>`,
:term:`fill patterns <fill pattern>`) behind calls such as
``fillMagneticGhosts(B, level, time)``.

The solver never talks to SAMRAI directly. When it needs ghost values, or
moments summed across patch borders, it asks the messenger of its level.


The picture
-----------

.. mermaid::

   flowchart TD
     MI["MultiPhysicsIntegrator<br/>one messenger per model pair"] --> IM["IMessenger"]
     IM --> HM["HybridMessenger"]
     IM --> MM["MHDMessenger<br/>MHD fine level, MHD coarse level"]
     HM --> HH["HybridHybridMessengerStrategy<br/>hybrid fine level, hybrid coarse level"]
     HM --> MH["MHDHybridMessengerStrategy<br/>hybrid fine level, MHD coarse level"]
     HH --> RP["refiners and synchronizers<br/>(what to move)"]
     RP --> SC["SAMRAI schedules, one per level<br/>(who sends what to whom)"]
     click HH href "messengers.html#the-hybrid-messenger" "The hybrid messenger"
     click SC href "messengers.html#when-schedules-are-built" "When schedules are built"

The messenger used for level ``L`` is chosen by the models of levels ``L-1``
and ``L``: its name is ``"HybridModel-HybridModel"``, ``"MHDModel-MHDModel"``
or ``"MHDModel-HybridModel"``. There is one messenger object per pair of
models, shared by all the levels that use that pair; level 0 uses the pair
(its model, its model) (``MultiPhysicsIntegrator::registerMessengers_``,
:src:`src/amr/multiphysics_integrator.hpp#L778-L819`). The C++ ``Simulator``
creates them once, at setup, through a ``MessengerFactory``
(``registerAndSetupMessengers``, called at
:src:`src/simulator/simulator.hpp#L324`).


How it works
------------

SAMRAI moves data in two layers. PHARE wraps both.

.. mermaid::

   flowchart TD
     A["algorithm: what to move<br/>built once, at setup"] --> S["schedule: who sends what to whom<br/>built per level, kept while levels last"]
     S --> F["fill: do the transfers<br/>every time the solver asks"]

#. **Algorithm: what to move.** A SAMRAI ``RefineAlgorithm`` (or
   ``CoarsenAlgorithm``) lists the quantities to move, with the operator that
   computes fine values from coarse ones and the fill pattern
   (:doc:`/dev/fill_patterns`). It does not depend on the patches. PHARE builds
   its algorithms once, at setup (``registerQuantities``,
   :src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L159-L203`).

#. **Schedule: who sends what to whom.** A ``RefineSchedule`` applies an
   algorithm to a given level: it knows every pair of patches that exchange
   data, on every MPI :term:`rank`. Building one is expensive and is
   *collective*: all ranks must build it together. It stays valid as long as the levels it
   was built from do not change.

#. **Fill.** ``schedule->fillData(time)`` does the transfers. This is what
   the messenger runs during the time step.

PHARE groups these in small classes (:src:`src/amr/messengers/refiner.hpp`,
:src:`src/amr/messengers/synchronizer.hpp`): a ``Refiner`` holds the
algorithms of one kind of transfer and one schedule per level; a
``RefinerPool`` holds the refiners of related quantities. Synchronizers are
the same for coarsening (fine to coarse, see :doc:`/dev/coarsening`).


The kinds of transfer
---------------------

Each ``Refiner`` has a type, ``RefinerType``
(:src:`src/amr/messengers/refiner.hpp#L16-L29`), which decides how its
schedule is built (``Refiner::makeSchedule_``,
:src:`src/amr/messengers/refiner.hpp#L60-L205`):

.. list-table::
   :header-rows: 1
   :widths: 26 44 30

   * - Type
     - What the schedule does
     - Used by (hybrid messenger)
   * - ``GhostField``
     - fills all ghosts: copies from the same level, refines from the coarser
       level where no same-level patch is (:doc:`/dev/ghost_cells`)
     - ``magGhostsRefiners_``, ``elecGhostsRefiners_`` (B, E)
   * - ``PatchGhostField``
     - copies from the same level only
     - ``chargeDensityPatchGhostsRefiners_``
   * - ``InitField``
     - fills a new level (cells and ghosts) by refining the coarser level
       everywhere, without reading the same level, which has no data yet
     - ``electricInitRefiners_`` (E on a new level)
   * - ``InitInteriorPart``
     - puts split coarse particles into the cells of a new level (not in its
       ghosts)
     - ``domainParticlesRefiners_``
   * - ``LevelBorderParticles``
     - puts split coarse particles in the layer just outside the level:
       level ghost particles (:doc:`/dev/level_ghost_particles`)
     - ``lvlGhostPartOldRefiners_``, ``lvlGhostPartNewRefiners_``
   * - ``ExteriorGhostParticles``
     - same-level particle exchange (:doc:`/dev/grabbing_particles`)
     - ``domainGhostPartRefiners_``
   * - ``PatchFieldBorderSum`` (and ``Vec``, ``Tensor``)
     - adds the neighbours' values on nodes in both patches' ghost boxes
       (:doc:`/dev/moments`)
     - ``popDensityBorderSumRefiners_``, ``popFluxBorderSumRefiners_``
   * - ``PatchFieldBorderMax`` (and ``Vec``)
     - same nodes, keeps the larger value
     - ``ionDensityBorderMaxRefiners_``, ``ionFluxBorderMaxRefiners_``
   * - ``LevelBorderField``
     - fills level ghosts only; not used at this commit
     - (none)


The hybrid messenger
--------------------

``HybridHybridMessengerStrategy``
(:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp`) is used
between two hybrid levels. Its calls during a time step:

.. list-table::
   :header-rows: 1
   :widths: 32 68

   * - Call
     - What it moves
   * - ``prepareStep``
     - nothing: saves B and E at the start of the step, for the finer level's
       time interpolation
   * - ``firstStep``
     - new level ghost particles, split from the coarser level
   * - ``fillMagneticGhosts``, ``fillElectricGhosts``
     - B and E ghosts (patch and level ghosts, interpolated in time)
   * - ``fillFluxBorders``, ``fillDensityBorders``
     - per-population moments summed across patch borders
   * - ``fillIonPopMomentGhosts``
     - no transfer: adds the level ghost particles' moments, old and new
       weighted by time
   * - ``fillIonGhostParticles``
     - particles that crossed into a neighbour patch
   * - ``fillIonBorders``
     - total charge density and flux made equal on shared border nodes
   * - ``lastStep``
     - no transfer: new level ghost particles become the old ones
   * - ``synchronize``, ``reflux``, ``postSynchronize``
     - fine to coarse at synchronization
       (:ref:`dev-one-time-step-sync`)

When each one is called is on :doc:`/dev/one_time_step`.


.. _dev-messengers-schedules:

When schedules are built
------------------------

.. mermaid::

   flowchart TD
     N["SAMRAI creates or replaces level L<br/>(init, regrid, rebalance)"] --> RL["messenger.registerLevel(L)<br/>same-level schedules of L"]
     RL --> INIT["level initializer fills L<br/>with init schedules, built and dropped"]
     INIT --> RH["SAMRAI: resetHierarchyConfiguration"]
     RH --> OH["messenger.onHierarchyChanged(L)<br/>for every level"]
     OH --> Q{"schedules of L built<br/>on the current levels?"}
     Q -->|yes| K["keep them"]
     Q -->|no| B["rebuild them"]

There are three lifetimes (``registerLevel`` and ``onHierarchyChanged``,
:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L208-L243`):

* **Same-level schedules** (patch ghosts, moment sums, particle exchange)
  read only level ``L``. They are built when ``L`` is created
  (``buildSameLevelSchedules_``,
  :src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L667-L700`).
* **Inter-level schedules** (B and E ghosts, level ghost particles,
  coarsening, refluxing) read ``L`` and the coarser levels. They are built in
  ``onHierarchyChanged``, after SAMRAI has finished replacing levels, and
  rebuilt only if one of the levels 0 to ``L`` is a different object from the
  one they were built on (``buildInterLevelSchedules_``,
  :src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L703-L722`).
* **Init-only schedules** (filling a new level) are used once: they are built
  at fill time and dropped (``Refiner::fillTransient``,
  :src:`src/amr/messengers/refiner.hpp#L219-L229`).

The calls come from ``MultiPhysicsIntegrator::initializeLevelData``
(:src:`src/amr/multiphysics_integrator.hpp#L327-L330`) and
``resetHierarchyConfiguration``
(:src:`src/amr/multiphysics_integrator.hpp#L347-L378`).


The tricks
----------

Schedules are rebuilt exactly when their levels change
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** An inter-level schedule of level ``L`` points at the coarser
levels it was built on. When level 0 is load-balanced, SAMRAI replaces level
0 and calls ``resetHierarchyConfiguration`` for level 0 only. A schedule of
level 2 that reads level 0 (through level 1, by recursive coarse
interpolation) would then point to a level that no longer exists. Rebuilding
everything at every regrid would be safe but costly.

**Fix.** Each messenger records, for each schedule set, the level objects it
was built from (``levelsOf_``, as weak pointers). ``onHierarchyChanged`` is
called for every level, whatever range SAMRAI passes, and rebuilds only the
schedules whose recorded levels are stale (``builtOnCurrentLevels_``,
:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L725-L744`; the
loop over all levels is at
:src:`src/amr/multiphysics_integrator.hpp#L372-L378`).

**History.** :pr:`1419` (inter-level schedules once per regrid; fixed a crash on
level-0 rebalance with 3 or more hybrid levels).

Level 0 has no coarser level, so its "ghost" refiners are same-level
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** The B and E ghost refiners are inter-level by nature, but level
0 also needs patch ghosts.

**Fix.** On the root level they are registered with the same-level schedules
(:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L680-L684`):
their schedule then only copies from the same level and applies the boundary
conditions.

Init schedules read the coarser level only
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** A new level has no data. A schedule that copies from the same
level first would copy garbage from the neighbour patches into the ghosts.

**Fix.** ``InitField`` schedules are built with no same-level source
(``nullptr``), so every value is refined from the coarser level
(:src:`src/amr/messengers/refiner.hpp#L144-L153`). Particles inside a new
level use ``PatchLevelInteriorFillPattern``, so that particle ghosts are never
filled by splitting: split particles placed independently in two patches'
ghost layers would not be exact copies of each other
(:src:`src/amr/messengers/refiner.hpp#L157-L170`).

After init and regrid: ghost times
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** Init and regrid schedules interpolate in space only, and read
the coarser levels at their latest time. A level created partway through a
coarser step gets ghost values from coarser data that is ahead of its own
time. For example: level 1 is at :math:`t + \Delta t/2`, level 0 already at
:math:`t + \Delta t`; a level 2 created now from level 1 may read level-0
values from :math:`t + \Delta t` in its ghosts.

**Fix (MHD).** The MHD messenger refills the ghosts of such a level with time
interpolation, as in a normal step, once its schedules are built
(``fillGhostsAtInitTime_`` and ``refillGhostsAtInitTime_``,
:src:`src/amr/messengers/mhd_messenger.hpp#L467-L494`). The hybrid messenger
has no such refill at this commit; its B and E ghosts are next filled by the
solver during the first sub-step.

**History.** :pr:`1494` (MHD ghost consistency).


Model pairs not yet coupled
---------------------------

``MHDHybridMessengerStrategy`` (a hybrid fine level over an MHD coarse level) exists so
that a hierarchy can mix models, but at this commit all its transfer methods
are empty (:src:`src/amr/messengers/mhd_hybrid_messenger_strategy.hpp`).

.. admonition:: In preparation
   :class: in-preparation

   Coupling a hybrid fine level to an MHD coarse level.


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
