.. _dev-one-time-step:

============================
One time step, top to bottom
============================

.. container:: breadcrumb

   You are here: :doc:`Developer guide </dev/index>` › One time step

**In one sentence:** one call to ``Simulator.advance()`` moves the whole
:term:`AMR` hierarchy forward by one coarse step :math:`\Delta t`: level 0
first, then each finer level in sub-steps, then the levels are made to agree.

This page is the top of the drill-down: every box links to the page that
explains it.


The call chain
--------------

.. flow::

   Python Simulator.advance() | then dump(): diagnostics, restarts | /dev/python_to_cpp
   calls
   C++ Simulator::advance() | chooses dt
   calls
   amr::Integrator::advance(dt) | may rebalance level 0 | /dev/load_balancing
   calls
   SAMRAI TimeRefinementIntegrator | steps the levels, coarsest first

   MultiPhysicsIntegrator::advanceLevel | once per step of each level
   calls
   messenger + solver (PPC or MHD) | ghosts, then the physics | /dev/messengers
   when levels meet
   synchronization | fine to coarse | /dev/coarsening
   every 2 steps of a level
   regrid finer levels | new patches where needed | /dev/tagging_regridding

* Python: ``advance()`` calls C++, then ``dump()`` writes restarts and
  diagnostics if due (:src:`pyphare/pyphare/simulator/simulator.py#L150-L165`,
  :src:`pyphare/pyphare/simulator/simulator.py#L230-L238`).
* ``Simulator::advance`` uses a fixed ``dt``, or an adaptive one (MHD only:
  the hybrid solver refuses it, :src:`src/amr/solvers/solver_ppc.hpp#L330`)
  (:src:`src/simulator/simulator.hpp#L558-L605`).
* SAMRAI owns the loop over levels and calls back PHARE's
  ``MultiPhysicsIntegrator``, its *level strategy*
  (:src:`src/simulator/simulator.hpp#L302`).
* Synchronization: ``standardLevelSynchronization``
  (:src:`src/amr/multiphysics_integrator.hpp#L566-L608`); regrid interval 2 is
  SAMRAI's default for a ratio of 2.


Levels are stepped coarsest first
---------------------------------

.. plot::
   :include-source: false
   :alt: Order in which SAMRAI advances three levels during one coarse step.

   from phare_draw import *
   from matplotlib.patches import Rectangle

   order = []
   def recurse(level, t0, t1):           # SAMRAI: one step of a level, then the finer ones
       order.append((level, t0, t1))
       if level < 2:
           dt = (t1 - t0) / 4
           for k in range(4):
               recurse(level + 1, t0 + k * dt, t0 + (k + 1) * dt)
   recurse(0, 0.0, 16.0)

   fig, (ax,) = strip(1, ["one coarse step, three levels: order of the 21 level steps"])
   col = [C["L0"], C["L1"], C["L2"]]
   for i, (lvl, t0, t1) in enumerate(order, start=1):
       y = 2 - lvl
       ax.add_patch(Rectangle((t0, y * 1.5), t1 - t0, 1.1, fc=col[lvl], ec="white", lw=1.5))
       text(ax, (t0 + t1) / 2, y * 1.5 + 0.55, str(i), color="white", bold=True,
            fs=13 if lvl < 2 else 11, halo=False)
   for lvl in range(3):
       text(ax, -0.4, (2 - lvl) * 1.5 + 0.55, f"level {lvl}", ha="right")
   for k, lab in enumerate(["t", "t+dt/4", "t+dt/2", "t+3dt/4", "t+dt"]):
       text(ax, 4 * k, -0.45, lab)
   ax.plot([-3.2, 16.5], [-0.8, 4.2], alpha=0)
   finish(fig, [ax])

Look at the numbers: level 0 first, then each level-1 step is followed by four
level-2 steps.

* A level twice finer takes 4 sub-steps: ``getMaxFinerLevelDt`` returns
  :math:`\Delta t / \text{ratio}^2` (:src:`src/amr/multiphysics_integrator.hpp#L420-L428`),
  for hybrid and MHD levels alike (ratio fixed at 2, :src:`src/amr/amr_constants.hpp#L9`).
* So when a fine level steps, its coarser level is **already at its new
  time**: fine ghosts are interpolated in time (next section).


Inside one level step
---------------------

.. plot::
   :include-source: false
   :alt: Coarse level saved at t and already at t+dt; a fine ghost at t_f is the time interpolation of the two.

   from phare_draw import *

   fig, (ax,) = strip(1, ["a fine ghost value between two coarse snapshots"])
   ax.plot([0, 8], [1.6, 1.6], color=C["L0"], lw=2)
   ax.plot([0, 8], [0, 0], color=C["L1"], lw=2)
   for k in range(5):
       ax.plot([2 * k, 2 * k], [-0.12, 0.12], color=C["L1"], lw=2)
   mark(ax, "primal", 0, 1.6, color=C["L0"], size=12)
   mark(ax, "primal", 8, 1.6, color=C["L0"], size=12)
   text(ax, 0, 2.15, "saved at t\n(prepareStep)")
   text(ax, 8, 2.15, "new value at t+dt\n(already computed)")
   mark(ax, "primal", 3, 0, color=C["B"], size=11)
   arrow(ax, (0.2, 1.45), (2.9, 0.15), color=C["grey"], lw=1.5)
   arrow(ax, (7.8, 1.45), (3.1, 0.15), color=C["grey"], lw=1.5)
   text(ax, 1.0, 0.75, r"$1-\alpha$", color=C["grey"])
   text(ax, 5.8, 0.95, r"$\alpha$", color=C["grey"])
   text(ax, 3, -0.6, r"fine ghost at $t_f$:  $\alpha = (t_f - t)/\Delta t = 0.375$", color=C["B"])
   text(ax, -0.7, 1.6, "level 0", ha="right")
   text(ax, -0.7, 0, "level 1", ha="right")
   ax.plot([-2.2, 9.6], [-0.9, 2.7], alpha=0)
   finish(fig, [ax])

Look at the fine ghost: it is a mix of the coarse value saved at :math:`t` and
the one already at :math:`t+\Delta t`. Here :math:`t_f` is the middle of the
second fine sub-step, where :math:`E_\text{avg}` is filled.

``MultiPhysicsIntegrator::advanceLevel``
(:src:`src/amr/multiphysics_integrator.hpp#L500-L561`), every time SAMRAI calls it:

* first sub-step only: ``messenger.firstStep`` makes the new level ghost
  particles (:doc:`/dev/level_ghost_particles`), ``solver.resetFluxSum``;
* ``solver.prepareStep`` and ``messenger.prepareStep`` save B (and E) at
  :math:`t` (:src:`src/amr/multiphysics_integrator.hpp#L533-L534`);
* ``solver.advanceLevel``: the physics (hybrid: :doc:`/dev/ppc_order`;
  MHD: :doc:`/dev/mhd_solver`);
* last sub-step only: ``messenger.lastStep``; finest level: accumulate the flux
  sum (not on level 0) and optional fine dump; then the load estimate.

The messenger is chosen by the models of the level and the coarser one; one
object per model pair (``registerMessengers_``,
:src:`src/amr/multiphysics_integrator.hpp#L778-L819`). Level 0, the *root*
level, also has one, for its patches only.


The tricks
----------

.. _dev-one-time-step-coarse-first:

Coarse level runs ahead, ghosts are interpolated in time
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** Fine ghosts are needed at times between two coarse steps.

**Fix.** ``prepareStep`` saves B and E at :math:`t`
(:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L580-L597`);
fine fills interpolate between them and the new coarse values
(:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L790-L831`), at
the time the data stands for (:math:`E_\text{avg}`: mid-sub-step,
:src:`src/amr/solvers/solver_ppc.hpp#L534-L562`). :pr:`1278`

Two saved copies of B
^^^^^^^^^^^^^^^^^^^^^

**Problem.** B at :math:`t` serves the finer level's interpolation and the
coarse level's refluxing.

**Fix.** The messenger and the solver each keep a copy
(:src:`src/amr/solvers/solver_ppc.hpp#L256-L272`), filled one after the other
(:src:`src/amr/multiphysics_integrator.hpp#L533-L534`). Do not merge them.

.. _dev-one-time-step-sync:

Coarse B is recomputed, not averaged (refluxing)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** The coarse B must follow the fine solution and keep div B = 0.

**Fix.** The solver sums :math:`E_\text{avg}/4` over the fine sub-steps
(:src:`src/amr/solvers/solver_ppc.hpp#L276-L288`; finest level in
:src:`src/amr/multiphysics_integrator.hpp#L549-L556`, others at sync
:src:`src/amr/multiphysics_integrator.hpp#L591-L600`); at sync this average
replaces the coarse :math:`E_\text{avg}` and coarse Faraday is redone
(:src:`src/amr/solvers/solver_ppc.hpp#L308-L320`). See :doc:`/dev/coarsening`.

Diagnostics come from Python, not from the level loop
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** ``MultiPhysicsIntegrator::dump_`` looks like the diagnostics
writer (:src:`src/amr/multiphysics_integrator.hpp#L464-L479`).

**Fix.** Regular dumps happen in ``Simulator.dump()`` after ``advance``
returns (:src:`pyphare/pyphare/simulator/simulator.py#L163`); ``dump_`` only
writes the optional fine dumps (``fine_dump_lvl_max``).


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
