.. _dev-subcycling:

============================================
Sub-cycling: finer levels take smaller steps
============================================

.. container:: breadcrumb

   You are here: :doc:`Developer guide </dev/index>` › :doc:`One time step </dev/one_time_step>` › Sub-cycling

**Applies to:** hybrid ✓ (fixed time step) | MHD ✓ (fixed or adaptive time step)

**In one sentence:** a refined level advances with several smaller time steps
for each step of the level below (*sub-cycling*), and after them the levels are
brought back to the same time and made to agree.

Why
---

A fine mesh needs a small time step. Hybrid whistler waves need
:math:`\Delta t \propto \Delta x^2`: halve the cell size and :math:`\Delta t`
must be divided by 4. With a fixed time step PHARE asks for exactly that bound;
with an adaptive one (MHD) each level asks for its own stable step. The
integrator (SAMRAI) then fills the coarse step with equal sub-steps.

.. claim C-sc-02

* ``getMaxFinerLevelDt`` returns the coarse step divided by ``ratio²``; the
  comment in the code gives the whistler reason
  (:src:`src/amr/multiphysics_integrator.hpp#L420-L428`). The ratio is 2
  (:src:`src/amr/amr_constants.hpp#L9`).

How many sub-steps?
-------------------

PHARE gives a bound on the next sub-step of a level. SAMRAI chooses the
number of sub-steps :math:`n`.

.. claim C-sc-03

* After each sub-step, SAMRAI takes the bound that PHARE's ``advanceLevel``
  returns, caps it at the previous bound times ``grow_dt`` (PHARE sets 1000, so
  that the sub-step can grow back), and sets :math:`n = \lceil \text{time
  remaining} / \text{bound} \rceil`, rounded up to a multiple of the regrid
  interval for a level that can be refined. The sub-step is *time remaining* /
  :math:`n`: the sub-steps left in a coarse step are equal
  (``TimeRefinementIntegrator::findNextDtAndStepsRemaining`` and
  ``advanceRecursivelyForRefinedTimestepping`` in SAMRAI, ``algs/
  TimeRefinementIntegrator.cpp``; SAMRAI at the commit pinned by the
  ``subprojects/samrai`` submodule, ``d246c078b``; ``grow_dt`` in the
  ``Integrator`` constructor, :pr:`1570`).

.. claim C-sc-01

* **Fixed time step (hybrid and MHD):** the bound is the level's fixed step,
  the coarse step divided by :math:`\text{ratio}^2`, so :math:`n = 4` at every
  level, sub-step :math:`\Delta t/4`. Hybrid only has a fixed step (the hybrid
  solver refuses an adaptive one, :src:`src/amr/solvers/solver_ppc.hpp#L330`).

.. claim C-sc-04

* **Adaptive time step (MHD):** the bound is the smaller of the coarse
  sub-step and the level's own stable step, computed on the level's new state
  after each sub-step (``MultiPhysicsIntegrator::advanceLevel`` and
  ``levelStableDt_``, :pr:`1570`). :math:`n` follows the fine level's own limit
  and can change from one coarse step to the next. Where the fine level sees the
  same wave speeds as the coarse one, :math:`n` is about 2 in ideal MHD, whose
  limit scales like :math:`\Delta x`, and up to 4 with Hall, whose whistler
  limit tends to :math:`\Delta x^2` on fine grids; faster waves on the fine
  level give more.
* The first sub-step of a level, after start-up or after the regrid that
  created it, has no stable step yet and uses the fixed-step bound, coarse
  step / :math:`\text{ratio}^2`. A uniform ideal-MHD state therefore takes 3
  sub-steps in the first coarse step (¼, then two of ⅜), then 2.
* Level 0 takes its own stable step (``computeStableDt``); the finer levels do
  not shorten it.
* The run stops with an error when a level's stable step is more than
  ``fine_level_substeps_max`` (default 1000) times smaller than its coarser
  level's, when the level-0 stable step is more than ``level0_dt_shrink_max``
  (default 1000) times smaller than the first one of the run (start or
  restart), or when a cell has a non-finite wave speed, the first of which the
  solver logs (``levelStableDt_``, :pr:`1570`).

Order of the advances
---------------------

.. plot::
   :include-source: false
   :alt: Order in which three levels advance during one coarse step, with the synchronizations.

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
   # synchronizations: level 2 -> level 1 at each end of a level-1 step, level 1 -> level 0 at the end
   for k in range(1, 5):
       x = 4 * k
       ax.plot([x, x], [1.12, 1.5], color=C["ok"], lw=3)
   ax.plot([16, 16], [2.62, 3.0], color=C["ok"], lw=3)
   text(ax, 16.4, 1.3, "2 → 1", ha="left", color=C["ok"])
   text(ax, 16.4, 2.8, "1 → 0", ha="left", color=C["ok"])
   for k, lab in enumerate([r"$t$", r"$t+\Delta t/4$", r"$t+\Delta t/2$", r"$t+3\Delta t/4$", r"$t+\Delta t$"]):
       text(ax, 4 * k, -0.45, lab)
   ax.plot([-3.2, 19.5], [-0.8, 4.2], alpha=0)
   finish(fig, [ax])

Look at the numbers, then at the green marks. Level 0 steps once; then each
level-1 step is followed by its four level-2 steps. With :math:`n = 4`
everywhere that is :math:`1 + 4 + 16 = 21` level steps.

.. claim C-sc-05

* SAMRAI advances a level, then recurses into the next finer level for *all* its
  sub-steps, then synchronizes (``TimeRefinementIntegrator.cpp``,
  ``advanceRecursivelyForRefinedTimestepping``).

Where the levels meet
---------------------

Each fine level is synchronized with the level below once per meeting time:
after every step of its coarser level. Level 2 meets level 1 after each of
level 1's four sub-steps; at the end of the coarse step level 1 also meets
level 0, finest pair first.

.. claim C-sc-06

* PHARE's ``standardLevelSynchronization`` loops from the finest level down:
  coarsen, reflux, refill ghosts (:src:`src/amr/multiphysics_integrator.hpp#L566-L608`).
  What each of these does is on :doc:`coarsening` and :doc:`refluxing`.
* At a *regrid point*, every second step of a level that can be refined, the
  same synchronizations happen, then the finer levels are regridded
  (:doc:`tagging_regridding`).

.. claim C-sc-07

Inside one level step
---------------------

.. flow::

   firstStep | first sub-step only: new level-ghost particles (hybrid); flux sums set to zero
   then
   prepareStep | save B, E (or the MHD state) at the start time | /dev/time_interpolation
   then
   solver.advanceLevel | the physics | /dev/ppc_order
   then
   lastStep | last sub-step only: level-ghost particles New → Old (hybrid)
   then
   accumulate flux sum | finest level only, weight sub-step / coarse step | /dev/refluxing

``MultiPhysicsIntegrator::advanceLevel``
(:src:`src/amr/multiphysics_integrator.hpp#L500-L562`) does these, in this
order, every time SAMRAI calls it.

.. claim C-sc-08

* ``firstStep``: hybrid fills the *New* level-ghost particles from the coarser
  level (not on level 0); MHD does nothing
  (:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L508-L530`,
  :src:`src/amr/messengers/mhd_messenger.hpp#L315-L321`).

.. claim C-sc-09

* The finest level adds its flux sums after each sub-step; every other level
  adds once, during the synchronization, after its own reflux
  (:src:`src/amr/multiphysics_integrator.hpp#L549-L556`,
  :src:`src/amr/multiphysics_integrator.hpp#L591-L600`).

.. claim C-sc-10

Next: :doc:`field_refinement`, how a coarse value becomes a fine one in space.

Go deeper
---------

* :doc:`time_interpolation`: what a fine ghost cell gets between two coarse times.
* :doc:`refluxing`: how the coarse level takes the fine answer at a meeting time.
* :doc:`tagging_regridding`: when and how the finer levels are rebuilt.
* :doc:`level_ghost_particles`: the particle side of the first and last sub-step.

Related
-------

* :doc:`one_time_step`: the whole call chain.
* :doc:`/numerics/amr`: the AMR algorithm for users.
