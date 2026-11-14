.. _dev-refluxing:

=================================================
Refluxing: the coarse level takes the fine answer
=================================================

.. container:: breadcrumb

   You are here: :doc:`Developer guide </dev/index>` › :doc:`Interpolation in time </dev/time_interpolation>` › Refluxing

**Applies to:** hybrid ✓ (the object that is summed is the electric field) | MHD ✓ (the objects are the fluxes)

**In one sentence:** at the time the levels meet, the coarse level repeats its
own update, this time using the electric field (hybrid) or the fluxes (MHD) that
the fine level computed along the coarse–fine border.

Why
---

The coarse level stepped first, across the coarse–fine border, with coarse
data. The fine level, with four times smaller steps, knows what really crossed
that border. Refluxing hands the coarse level this better answer, so that the
coarse level agrees with the fine level at the border (conservation, see the
guarantees below), and :math:`\nabla\cdot B` stays zero.

Three phases
------------

1. **Accumulate.** After each fine sub-step the fine level adds a quarter of its
   own flux (hybrid: of :math:`E_\text{avg}`) to a running sum.
2. **Coarsen.** At the meeting time the sum is copied or averaged onto the
   coarse faces (or edges) it covers.
3. **Redo.** The coarse level repeats its update over the whole coarse step,
   from the values it saved at the start of the step, using the new fluxes.

.. claim C-rx-02

* Each level's sums are zeroed at its first sub-step. The finest level adds
  :math:`1/\text{ratio}^2 = 1/4` times its quantity after every sub-step
  (hybrid: :math:`E_\text{avg}`, MHD: the step's fluxes and E)
  (:src:`src/amr/multiphysics_integrator.hpp#L526-L532`,
  :src:`src/amr/multiphysics_integrator.hpp#L550-L556`,
  :src:`src/amr/solvers/solver_ppc.hpp#L277-L298`,
  :src:`src/amr/solvers/solver_mhd.hpp#L310-L342`).
* A level that has finer levels adds, once, after its own reflux, so that the
  contribution of all finer levels travels up the chain
  (:src:`src/amr/multiphysics_integrator.hpp#L591-L600`).

.. claim C-rx-04

When it runs
------------

.. flow::

   synchronize | hybrid: coarsen E, charge density, charge flux. MHD: nothing | /dev/coarsening
   then
   messenger.reflux | coarsen the summed flux onto the coarse level, refill coarse patch ghosts
   then
   solver.reflux | the coarse level redoes its update
   then
   accumulate | for the next level up, if the coarse level is not level 0
   then
   postSynchronize | refill patch ghosts (hybrid)

This is one pass of ``standardLevelSynchronization``, finest pair first
(:src:`src/amr/multiphysics_integrator.hpp#L575-L609`); when the levels meet is
on :doc:`subcycling`.

.. claim C-rx-01

Hybrid
------

.. plot::
   :include-source: false
   :alt: One dimension. Four fine values of the electric field at the border node are averaged into one coarse value, which the coarse Faraday update then uses for coarse cell 1.

   from phare_draw import *
   from matplotlib.patches import Rectangle

   fig, (ax,) = strip(1, ["hybrid, 1D: the border node x = 2"], aspect=False, height=3.4)
   # coarse cells 0..4 (dx = 1), fine cells 4..11 under coarse 2..5 (dx = 1/2)
   for c in range(5):
       ax.add_patch(Rectangle((c, 2.2), 1, 0.6, fc=C["faint"], ec=C["L0"], lw=1.8))
   for f in range(4, 10):
       ax.add_patch(Rectangle((f / 2, 0.5), 0.5, 0.6, fc=C["patch"], ec=C["L1"], lw=1.8))
   text(ax, -0.15, 2.5, "level 0", ha="right", color=C["L0"])
   text(ax, -0.15, 0.8, "level 1", ha="right", color=C["L1"])
   # the border node x = 2: four fine sub-step values, one coarse value
   for k in range(4):
       ax.add_patch(Rectangle((2 - 0.06 + 0.03 * k, 1.25 + 0.17 * k), 0.12, 0.13, fc=C["E"], ec="white"))
   text(ax, 2.25, 1.55, r"$E_1 \ldots E_4$", ha="left", color=C["E"])
   mark(ax, "primal", 2, 2.5, color=C["E"], size=12)
   arrow(ax, (2, 1.9), (2, 2.3), color=C["E"], lw=2)
   text(ax, 2, 3.25, r"coarse $E_\mathrm{avg}$ at $x=2$ = $\frac{1}{4} \sum E_k$", color=C["E"])
   text(ax, 0.9, 1.7, "cell 1:\nits B changes", color=C["B"], ha="center")
   ax.plot([-1.2, 5.6], [-0.2, 3.7], alpha=0)
   finish(fig, [ax])

Look at the border node: the four sub-step values of the fine level become one
value on the coarse level; the coarse cell next to it, cell 1, has its magnetic
field redone with it.

1. **Accumulate:** each fine sub-step adds ¼·:math:`E_\text{avg}`, the field
   that its corrector Faraday step used
   (:src:`src/amr/solvers/solver_ppc.hpp#L277-L288`).

.. claim C-rf-06

2. **Coarsen:** the sum replaces the coarse :math:`E_\text{avg}` over the nodes
   or edges the fine level covers, with the same stencil as
   :doc:`coarsening` (average of 2 fine edges along a dual direction, copy
   of the coincident fine node otherwise), then the coarse patch ghosts are
   refilled
   (:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L189-L199`,
   :src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L641-L646`).

.. claim C-rx-05

3. **Redo:** the coarse solver computes :math:`B = B_\text{old} - \Delta t\,\nabla\times E_\text{avg}` over the
   *whole* coarse level, with :math:`B_\text{old}` its own copy of B saved at the
   start of the coarse step and :math:`\Delta t` the whole coarse step, then
   refills the B ghosts
   (:src:`src/amr/solvers/solver_ppc.hpp#L257-L272`,
   :src:`src/amr/solvers/solver_ppc.hpp#L309-L321`).

.. claim C-rx-06

* Away from the fine level the redo gives the same B as the coarse corrector
  did: the same starting B (the predictors write another array) and the same
  :math:`E_\text{avg}` (:src:`src/amr/solvers/solver_ppc.hpp#L390-L403`,
  :src:`src/amr/solvers/solver_ppc.hpp#L498-L503`).

.. claim C-rf-08

* The coarse B is never coarsened from the fine B: B is not in any
  synchronizer, only E, charge density and charge flux are
  (:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L912-L924`).
  The coarse edge E is an average of fine edges or a copy, so the change of the
  flux through a coarse face is the sum of the changes of its fine faces
  (:src:`src/amr/data/field/coarsening/electric_field_coarsener.hpp#L45-L145`).
  Separately, :math:`\nabla\cdot B` is unchanged whatever E is, because the
  divergence of a discrete curl is zero.

.. claim C-rf-09 C-rf-10 C-rx-08

MHD
---

.. plot::
   :include-source: false
   :alt: One dimension. Four sub-step fluxes at the border face are averaged into one flux that replaces the coarse flux there.

   from phare_draw import *
   from matplotlib.patches import Rectangle

   fig, (ax,) = strip(1, ["MHD, 1D, fixed time step: the border face x = 2"], aspect=False, height=3.4)
   for c in range(5):
       ax.add_patch(Rectangle((c, 2.2), 1, 0.6, fc=C["faint"], ec=C["L0"], lw=1.8))
   for f in range(4, 10):
       ax.add_patch(Rectangle((f / 2, 0.5), 0.5, 0.6, fc=C["patch"], ec=C["L1"], lw=1.8))
   text(ax, -0.15, 2.5, "level 0", ha="right", color=C["L0"])
   text(ax, -0.15, 0.8, "level 1", ha="right", color=C["L1"])
   for k in range(4):
       ax.add_patch(Rectangle((2 - 0.06 + 0.03 * k, 1.25 + 0.17 * k), 0.12, 0.13, fc=C["E"], ec="white"))
   text(ax, 2.25, 1.55, r"$F_1 \ldots F_4$", ha="left", color=C["E"])
   arrow(ax, (2, 1.9), (2, 2.3), color=C["E"], lw=2)
   text(ax, 2, 3.25, r"flux through $x=2$ = $\frac{1}{4} \sum F_k$", color=C["E"])
   text(ax, 0.9, 1.7, "cell 1 gets\nwhat the fine\ncells receive", color=C["B"], ha="center")
   ax.plot([-1.2, 5.6], [-0.2, 3.7], alpha=0)
   finish(fig, [ax])

Look at the border face: with a fixed time step, the mass leaving coarse cell 1
through :math:`x = 2` is the mass the fine cells received.

1. **Accumulate:** each fine sub-step adds ¼·its step flux, a time average over
   that sub-step (the Runge–Kutta stage fluxes are weighted, the weights sum to 1 for all four integrators)
   (:src:`src/amr/solvers/solver_mhd.hpp#L310-L327`,
   :src:`src/amr/solvers/time_integrator/tvdrk2_integrator.hpp#L35-L40`).

.. claim C-rx-02 C-rf-01

2. **Coarsen:** the sums of the :math:`\rho`, :math:`\rho \mathbf v` and
   :math:`E_\text{tot}` fluxes are coarsened onto the coarse step fluxes (1D
   copy, 2D mean of 2 fine faces, 3D mean of 4), and the summed E onto the
   coarse E; then the patch ghosts of the fluxes are refilled
   (:src:`src/amr/messengers/mhd_messenger.hpp#L362-L367`,
   :src:`src/amr/messengers/mhd_messenger.hpp#L601-L615`,
   :src:`src/amr/data/field/coarsening/mhd_flux_coarsener.hpp#L46-L149`).

.. claim C-rx-05 C-rf-03

3. **Redo:** the coarse level repeats **one** Euler update of the whole level
   from the state it saved at the start of the step, with these fluxes and E
   (B through Faraday with that E), over :math:`\Delta t` = sync time − step start,
   then refills the state ghosts
   (:src:`src/amr/solvers/solver_mhd.hpp#L276-L306`,
   :src:`src/amr/solvers/solver_mhd.hpp#L345-L357`,
   :src:`src/amr/solvers/time_integrator/euler_using_computed_flux.hpp#L24-L35`).

.. claim C-rx-07

* MHD ``synchronize`` is empty: no fine state is averaged onto the coarse state
  (:src:`src/amr/messengers/mhd_messenger.hpp#L360`). The coarse cells under the
  fine level are *recomputed* by the redo, from the coarsened fine fluxes on
  every covered face; they equal the mean of their fine cells only because the
  fine level started from a refinement that keeps the coarse average
  (:doc:`field_refinement`).

.. claim C-rx-08 C-fr-12

What it guarantees
------------------

* Hybrid and MHD: the coarse level uses the fine answer at the border.
* MHD with a fixed time step: the flux the outside coarse cell sees equals the
  fine time-integrated flux, because four sub-steps of :math:`\Delta t/4`,
  each weighted by ¼, make the mean flux, and the mean times :math:`\Delta t`
  is the fine time integral.
* MHD, any time step: the level-0 total changes only through physical
  boundaries (each coarse face has one flux, used by both its cells).

.. claim C-rf-02 C-rf-05

.. admonition:: Known issue (MHD, adaptive time step)
   :class: note

   The weight ¼ per sub-step is fixed (:src:`src/amr/multiphysics_integrator.hpp#L553-L555`).
   With an adaptive step the number of sub-steps can differ from 4
   (:doc:`subcycling`), and the sum is then no longer the mean flux. The coarse
   and fine levels then no longer match at the border.

.. claim C-rx-03 C-sc-04

Next: :doc:`coarsening`, what else goes from fine to coarse at the meeting time.

Go deeper
---------

* :doc:`subcycling`: when the levels meet.
* :doc:`coarsening`: the coarsening stencils.
* :doc:`ppc_order`: the hybrid step, with the two Faraday updates.
* :doc:`mhd_solver`: the MHD step and its stage fluxes.

Related
-------

* :doc:`field_refinement`: how the fine level starts from the coarse one.
* :doc:`/numerics/amr`: the AMR algorithm for users.
