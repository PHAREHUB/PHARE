.. _dev-one-time-step:

============================
One time step, top to bottom
============================

.. container:: breadcrumb

   You are here: :doc:`Developer guide </dev/index>` › One time step

**In one sentence:** one call to ``Simulator.advance()`` in Python moves the
whole AMR (adaptive mesh refinement) hierarchy of levels forward by one coarse time step :math:`\Delta t`: the
coarsest level first, then each finer level in smaller sub-steps, then the
levels are made to agree again.

This page is the top of the drill-down. Every box in the diagrams below is a
link to the page that explains it.


The picture
-----------

.. mermaid::

   flowchart TD
     PY["Python<br/>Simulator.advance()"] --> CPP["C++ Simulator::advance()<br/>chooses dt"]
     PY -->|"then"| DIAG["dump(): diagnostics<br/>and restarts"]
     CPP --> INT["amr::Integrator::advance"]
     INT --> TRI["SAMRAI<br/>TimeRefinementIntegrator"]
     TRI -->|"each level step"| AL["MultiPhysicsIntegrator<br/>::advanceLevel"]
     TRI -->|"levels meet"| SYNC["Synchronization"]
     TRI -->|"every 2 steps"| RG["Regrid<br/>finer levels"]
     AL -->|"hybrid level"| PPC["SolverPPC"]
     AL -->|"or MHD level"| MHD["SolverMHD"]
     PPC --> MSG["Messenger"]
     MHD --> MSG
     click PY href "python_to_cpp.html" "From Python script to C++"
     click MSG href "messengers.html" "Messengers and schedules"
     click PPC href "ppc_order.html" "Push and deposition order in the PPC solver"
     click MHD href "mhd_solver.html" "MHD solver internals"
     click SYNC href "coarsening.html" "Coarsening"
     click RG href "tagging_regridding.html" "Tagging and regridding"
     click DIAG href "diagnostics_restarts.html" "Diagnostics and restarts internals"

Words used in the picture:

* **Messenger**: the object that moves data between a level and the next
  coarser one, and between neighbouring patches of the same level.
* **Ghost cells**: cells just outside a patch, filled with copies or
  interpolations of data owned by a neighbour or by the coarser level
  (:doc:`/dev/ghost_cells`).
* **Level ghost particles**: particles placed just outside a fine level,
  made from coarser-level particles (:doc:`/dev/level_ghost_particles`).
* **Synchronization**: once a fine level has caught up in time, its better
  solution is copied (coarsened) onto the coarser level underneath.
* **Reflux**: the part of synchronization that redoes the coarse magnetic
  field update with the fine electric field (see the tricks below).
* **Regrid**: rebuild a finer level with new patches where refinement is
  now needed.
* Each level runs one solver, chosen by the model of that level: the hybrid
  solver ``SolverPPC`` (PPC: predictor-predictor-corrector) or the MHD
  (magnetohydrodynamics, a fluid model) solver ``SolverMHD``.


How it works
------------

#. **Python.** ``Simulator.advance()``
   (:src:`pyphare/pyphare/simulator/simulator.py#L150-L165`) calls the C++
   simulator (``cpp_sim.advance()``, a C++ function exposed to Python with
   the pybind11 library), then calls ``dump()``,
   which writes restarts and diagnostics if one is due at the new time
   (:src:`pyphare/pyphare/simulator/simulator.py#L230-L238`).
   ``Simulator.run()`` is just a loop over ``advance()``.

#. **C++ simulator.** ``Simulator::advance``
   (:src:`src/simulator/simulator.hpp#L558-L605`) uses a fixed ``dt``, or,
   with an adaptive time step, asks every level for its stable ``dt``,
   converts it to the equivalent level-0 step (times 4 per level, since each
   finer level takes 4 sub-steps, see below), and keeps the smallest over all
   levels and MPI ranks (``MultiPhysicsIntegrator::computeStableDt``,
   :src:`src/amr/multiphysics_integrator.hpp#L435-L461`). Then it calls
   ``integrator_->advance(dt_)``.

#. **AMR integrator.** ``amr::Integrator::advance``
   (:src:`src/amr/wrappers/integrator.hpp#L57-L64`) decides whether the
   coarsest level must be load-balanced now, and hands over to SAMRAI's
   ``TimeRefinementIntegrator::advanceHierarchy``.

#. **SAMRAI steps the levels.** SAMRAI (the AMR library PHARE is built on)
   owns the loop over levels and sub-steps. It does not know any physics: for
   each step of each level it calls back PHARE's ``MultiPhysicsIntegrator``.
   That class is SAMRAI's *level strategy*: the object SAMRAI calls to
   advance, synchronize or initialize one level
   (:src:`src/simulator/simulator.hpp#L302`). The order of the calls is shown
   in the next picture.

#. **One level, one step.** ``MultiPhysicsIntegrator::advanceLevel``
   (:src:`src/amr/multiphysics_integrator.hpp#L500-L561`) picks the solver,
   the model and the messenger of that level and runs the steps of the
   "inside one level step" diagram below.

#. **Synchronization.** When a fine level has caught up with its coarser
   level, SAMRAI calls ``MultiPhysicsIntegrator::standardLevelSynchronization``
   (:src:`src/amr/multiphysics_integrator.hpp#L566-L608`): fine data is
   coarsened onto the coarser level, finest level first (see
   :ref:`dev-one-time-step-sync`).

#. **Regridding.** Every 2 steps of a level, SAMRAI may rebuild the finer
   levels (new patches where refinement is needed). The interval 2 is not
   set by PHARE: SAMRAI's ``TimeRefinementIntegrator`` derives it from the
   refinement ratio (2 for an even ratio), and only regrids levels that are
   allowed to have a finer level. See :doc:`/dev/tagging_regridding`.


Levels are stepped coarsest first
---------------------------------

The refinement ratio between two levels is always 2
(``refinementRatio``, :src:`src/amr/amr_constants.hpp#L9`). The hybrid
equations have whistler waves, whose stable time step scales like
:math:`\Delta x^2`. So a level twice finer needs a time step four times
smaller: ``MultiPhysicsIntegrator::getMaxFinerLevelDt`` returns
:math:`\Delta t_\text{coarse} / \text{ratio}^2 = \Delta t_\text{coarse} / 4`
(:src:`src/amr/multiphysics_integrator.hpp#L420-L428`). With a fixed time
step, each level therefore takes **4 sub-steps** per step of the level above.
This rule does not depend on the model: MHD levels are sub-stepped by 4 too.

SAMRAI walks the hierarchy depth first: it advances a level by one step, then
advances all finer levels up to the same time, then moves on.

.. plot::
   :include-source: false
   :alt: Order in which SAMRAI advances three levels during one coarse step.

   import matplotlib.pyplot as plt
   from matplotlib.patches import Rectangle

   n_levels, ratio_t = 3, 4   # 3 levels, 4 sub-steps per coarse step
   order = []

   def recurse(level, t0, t1):
       # SAMRAI: for each step of `level`: advance it, then the finer levels
       order.append((level, t0, t1))
       if level + 1 < n_levels:
           dt = (t1 - t0) / ratio_t
           for k in range(ratio_t):
               recurse(level + 1, t0 + k * dt, t0 + (k + 1) * dt)

   recurse(0, 0.0, 1.0)

   colors = ["#4C72B0", "#55A868", "#C44E52"]
   fig, ax = plt.subplots(figsize=(6.5, 3.4))
   for i, (lvl, t0, t1) in enumerate(order, start=1):
       y = n_levels - 1 - lvl
       ax.add_patch(Rectangle((t0, y + 0.15), t1 - t0, 0.7,
                              facecolor=colors[lvl], edgecolor="white", lw=1.5))
       ax.text((t0 + t1) / 2, y + 0.5, str(i), ha="center", va="center",
               color="white", fontsize=14 if lvl < 2 else 11, fontweight="bold")

   # synchronizations: L2 -> L1 at the end of each L1 step, everything at t + dt
   for k in range(1, 5):
       t = k / 4
       ax.annotate("", xy=(t, 1.15), xytext=(t, 0.85),
                   arrowprops=dict(arrowstyle="-|>", color="0.3", lw=1.5))
   ax.annotate("", xy=(1.0, 2.15), xytext=(1.0, 1.85),
               arrowprops=dict(arrowstyle="-|>", color="0.3", lw=1.5))
   ax.text(1.01, 1.0, "sync", va="center", fontsize=11, color="0.3")
   ax.text(1.01, 2.0, "sync", va="center", fontsize=11, color="0.3")

   ax.set_xlim(-0.02, 1.1)
   ax.tick_params(labelsize=11)
   ax.set_ylim(0, n_levels)
   ax.set_yticks([n_levels - 1 - l + 0.5 for l in range(n_levels)])
   ax.set_yticklabels([f"level {l}" for l in range(n_levels)])
   ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
   ax.set_xticklabels(["t", "t+dt/4", "t+dt/2", "t+3dt/4", "t+dt"])
   for s in ("top", "right", "left"):
       ax.spines[s].set_visible(False)
   ax.tick_params(axis="y", length=0)
   fig.tight_layout()

*One coarse step with three levels. Numbers give the order of the 21 calls to
``advanceLevel``. Level 0 jumps to* :math:`t+\Delta t` *first (call 1). Level
1 then does 4 sub-steps, and after each one, level 2 does 4 sub-steps to catch
up. Arrows mark the synchronizations: level 2 onto level 1 at the end of each
level-1 step, and all levels at* :math:`t+\Delta t`. *Regridding (not drawn)
also happens at some of these points.*

This order has one big consequence: **when a fine level is advanced, its
coarser level is already at its new time.** The fine level never waits for
data from the future. It interpolates in time between two coarse snapshots it
already has (see :ref:`dev-one-time-step-coarse-first`).


Inside one level step
---------------------

What ``MultiPhysicsIntegrator::advanceLevel``
(:src:`src/amr/multiphysics_integrator.hpp#L500-L561`) does every time SAMRAI
calls it. "Messenger" is the object that moves data between this level and the
next coarser one (and between patches of this level).

.. mermaid::

   flowchart TD
     A{"first sub-step of this level<br/>in this coarse step?"}
     A -->|yes| B["messenger.firstStep<br/>make the new level ghost particles<br/>from the coarser level, already at its new time"]
     B --> C["solver.resetFluxSum"]
     A -->|no| D
     C --> D["solver.prepareStep and messenger.prepareStep<br/>keep copies of B and E at time t"]
     D --> E["solver.advanceLevel<br/>t to t + dt on every patch of the level"]
     E --> F{"last sub-step?"}
     F -->|yes| G["messenger.lastStep<br/>new level ghost particles become the old ones"]
     F -->|no| H
     G --> H{"finest level?"}
     H -->|yes| I["optional fine-level dump<br/>accumulate the flux sum (not on level 0)"]
     H -->|no| J
     I --> J["load-balancer workload estimate"]
     click B href "level_ghost_particles.html" "Level ghost particles"
     click G href "level_ghost_particles.html" "Level ghost particles"
     click D href "messengers.html" "Messengers and schedules"
     click E href "ppc_order.html" "Push and deposition order in the PPC solver"

The messenger used depends on the models of the two levels: ``HybridModel-HybridModel``,
``MHDModel-MHDModel`` or ``MHDModel-HybridModel``. There is one messenger
object per pair of models, shared by all levels that use that pair
(``registerMessengers_``, :src:`src/amr/multiphysics_integrator.hpp#L778-L819`).
Level 0, also called the *root* level in the code, has a messenger too: there, it only moves data between patches of
level 0.

Inside the hybrid solver
^^^^^^^^^^^^^^^^^^^^^^^^

``SolverPPC::advanceLevel`` (:src:`src/amr/solvers/solver_ppc.hpp#L378-L404`)
advances one level by one step with a :term:`predictor-predictor-corrector
<PPC>` scheme. Particles
are pushed twice: the first push only *predicts* the ion moments at the new
time (particles stay where they were); the second push really moves them.

.. mermaid::

   flowchart TD
     P1["predictor 1<br/>Faraday, Ampère, Ohm<br/>gives predicted B, E"] --> A1["average<br/>(now + predicted) / 2"]
     A1 --> M1["move ions, moments only<br/>push a copy, deposit density, flux"]
     M1 --> P2["predictor 2<br/>Faraday, Ampère, Ohm<br/>with predicted moments"]
     P2 --> A2["average"]
     A2 --> M2["move ions, for real<br/>push, exchange, deposit"]
     M2 --> C["corrector<br/>Faraday, Ampère, Ohm<br/>gives B, E at the new time"]
     click M1 href "ppc_order.html" "Push and deposition order"
     click M2 href "ppc_order.html" "Push and deposition order"

Each "move ions" step (``SolverPPC::moveIons_``,
:src:`src/amr/solvers/solver_ppc.hpp#L590-L666`) does:

.. mermaid::

   flowchart TD
     U["push and deposit<br/>(IonUpdater)"] --> B["sum moments<br/>over patch borders"]
     B --> L["add level ghost<br/>particles' moments"]
     L --> G["second move only:<br/>grab particles from<br/>neighbour patches"]
     G --> I["total charge density,<br/>bulk velocity, charge flux"]
     I --> F["make these totals agree<br/>on patch borders"]
     click U href "ppc_order.html" "Push and deposition order"
     click B href "moments.html" "Moments: deposition and sum borders"
     click L href "level_ghost_particles.html" "Level ghost particles"
     click G href "grabbing_particles.html" "Grabbing particles from neighbours"
     click F href "moments.html" "Moments: deposition and sum borders"

The full order, and why it is this order, is on :doc:`/dev/ppc_order`.


.. _dev-one-time-step-sync:

Synchronization: making the levels agree
----------------------------------------

Once a fine level has reached the time of its coarser level, the coarse level
must "feel" the better fine solution where the two overlap.
``standardLevelSynchronization``
(:src:`src/amr/multiphysics_integrator.hpp#L566-L608`) goes from the finest
level down. For each fine level, with the hybrid messenger:

.. mermaid::

   flowchart TD
     S1["synchronize<br/>coarsen E, ion charge density<br/>and ion charge flux"] --> S2["messenger.reflux<br/>coarsen the time-averaged<br/>fine E onto the coarse level"]
     S2 --> S3["coarse solver.reflux<br/>redo coarse Faraday<br/>from B at the start of the step"]
     S3 --> S4["postSynchronize<br/>refill coarse patch ghosts<br/>of E and charge density"]
     click S1 href "coarsening.html" "Coarsening"
     click S2 href "coarsening.html" "Coarsening"
     click S4 href "ghost_cells.html" "Ghost cells"

The functions are ``HybridHybridMessengerStrategy::synchronize``, ``reflux``
and ``postSynchronize``
(:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L626-L663`) and
``SolverPPC::reflux`` (:src:`src/amr/solvers/solver_ppc.hpp#L308-L320`).


The tricks
----------

.. _dev-one-time-step-coarse-first:

The coarse level runs ahead; fine ghosts are interpolated in time
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** A fine level needs values in its ghost cells along the
coarse-fine boundary at every sub-step, at times between two coarse steps.

**Fix.** Before a level is advanced, ``HybridHybridMessengerStrategy::prepareStep``
(:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L580-L597`)
saves its B and E at the start time (``Bold_``, ``Eold_``). When the finer
level later fills its coarse-fine ghosts at a time :math:`t_f`, the values are
interpolated linearly in time between that saved copy and the coarse level's
new B and E (registration at
:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L790-L831`). The
solver asks for each fill at the time the data stands for: for example, the
average field :math:`E_\text{avg}` is filled at the middle of the fine
sub-step, :math:`(t + t_\text{new})/2`
(``SolverPPC::average_``, :src:`src/amr/solvers/solver_ppc.hpp#L534-L562`).
Particles use the same idea with two sets of level ghost particles, old and
new (:doc:`/dev/level_ghost_particles`).

**History.** :pr:`1278` (second-order coarse-fine time interpolation) made B time
interpolated too and set the fill times to the solver's sub-step times.

Two saved copies of B at the start of the step
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** The start-of-step magnetic field is needed for two different
reasons: the finer level interpolates its ghosts in time from it, and the
coarse level redoes its Faraday update from it at synchronization (refluxing).

**Fix.** Both keep their own copy. The messenger copies B and E in
``prepareStep`` (above), and the solver copies B in ``SolverPPC::prepareStep``
(:src:`src/amr/solvers/solver_ppc.hpp#L256-L272`). They are filled one after
the other in ``advanceLevel``
(:src:`src/amr/multiphysics_integrator.hpp#L533-L534`). Do not merge them
without checking both uses.

Coarse B is recomputed, not averaged from fine B (refluxing)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** Where a fine level covers a coarse one, the coarse magnetic field
should match the fine solution, and keep :math:`\nabla\cdot B = 0`.

**Fix.** The coarse B is not coarsened from fine B. During each fine sub-step,
the solver adds :math:`E_\text{avg}/4` to a running sum ``fluxSumE_``
(``SolverPPC::accumulateFluxSum``,
:src:`src/amr/solvers/solver_ppc.hpp#L276-L288`), so after 4 sub-steps it
holds the time-averaged fine electric field. At synchronization this average
is coarsened onto the coarse :math:`E_\text{avg}`, and the coarse level redoes
its Faraday update :math:`B = B_\text{old} - \Delta t\,\nabla\times E_\text{avg}`
(``SolverPPC::reflux``). B changes by a discrete curl, so its discrete
divergence does not change.

The finest level accumulates its sum after each sub-step, inside
``advanceLevel`` (:src:`src/amr/multiphysics_integrator.hpp#L549-L556`).
Intermediate levels can only do it at synchronization, after their own E
has been corrected by the finer level
(:src:`src/amr/multiphysics_integrator.hpp#L591-L600`).

Regular diagnostics are written from Python, not from the level loop
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** ``MultiPhysicsIntegrator::dump_``
(:src:`src/amr/multiphysics_integrator.hpp#L464-L479`) looks like the place
where diagnostics are written. It is not.

**Fix.** The regular diagnostics and restarts are written by
``Simulator.dump()`` in Python, after ``cpp_sim.advance()`` returns, so always
at a coarse time (:src:`pyphare/pyphare/simulator/simulator.py#L163`).
``dump_`` only writes the optional *fine dumps* (``fine_dump_lvl_max``):
fine-level data at sub-step times between two coarse steps.


Go deeper
---------

* :doc:`/dev/messengers`: what each messenger fills, and when its schedules are built.
* :doc:`/dev/ghost_cells`: patch ghost cells versus level ghost cells.
* :doc:`/dev/level_ghost_particles`: old and new level ghost particles.
* :doc:`/dev/ppc_order`: the hybrid solver step by step.
* :doc:`/dev/coarsening`: what synchronization coarsens, and how.
* :doc:`/dev/tagging_regridding`: when and how the finer levels are rebuilt.

Related
-------

* :doc:`/dev/python_to_cpp`: how the Python script becomes the C++ simulator that this page runs.
* :doc:`/dev/mhd_solver`: the MHD solver step.
* :doc:`/dev/diagnostics_restarts`: what ``dump()`` writes.
* :doc:`/numerics/amr`: the AMR algorithm for users.
