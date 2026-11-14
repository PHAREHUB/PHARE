================
Time integration
================

**Goal:** know how PHARE advances the hybrid and MHD equations over one time step, and
how to choose the time step.

.. contents:: On this page
   :local:
   :depth: 1

Setting the time step
---------------------

The time step is set in the ``Simulation`` block. These lines are a fragment of a full
input script (for a complete one, see :doc:`/getting_started/first_simulation`):

.. code-block:: python

   # fixed time step: give two of time_step, time_step_nbr, final_time
   sim = ph.Simulation(
       time_step=0.005,
       final_time=10.0,
       # ... other options
   )

   # MHD only: time step recomputed at every step from a stability condition
   sim = ph.Simulation(
       time_step={"mode": "adaptive", "cfl_wave": 0.8},
       final_time=10.0,
       # ... other options
   )

The time step you give is the one of the coarsest :term:`level`. Each finer level takes smaller
steps (see `Time steps on refined levels`_).

The hybrid scheme: predictor, predictor, corrector
--------------------------------------------------

(:term:`PPC` for short.)

The difficulty in a hybrid code is that the electric field depends on the ion moments at
the *end* of the step, which are only known after the particles have moved, and the
particles move in that same field. PHARE solves this with two predictions and one
correction. Starting from :math:`\mathbf{B}^n`, :math:`\mathbf{E}^n` and the ion moments
at time :math:`n`:

.. plot::
   :include-source: False
   :caption: Red: fields; green: particles and ion moments. Each stage starts again from Bⁿ and
      advances with a mid-step (averaged) field; only the second push is kept.

   from phare_draw import *

   T = ["Predictor 1, then a trial push",
        "Predictor 2, then the real push",
        "Corrector"]
   rows = [("B*, E*  (Ohm with nⁿ)", "ñⁿ⁺¹: predicted, particles discarded", "Faraday with Eⁿ"),
           ("B**, E**  (Ohm with ñⁿ⁺¹)", "nⁿ⁺¹: particles kept", "Faraday with average 1"),
           ("Bⁿ⁺¹, Eⁿ⁺¹  (Ohm with nⁿ⁺¹)", None, "Faraday with average 2")]
   fig, ax = strip(3, titles=T, aspect=False, height=1.6)
   for k, (a, (flab, push, far)) in enumerate(zip(ax, rows)):
       for x, lab in [(0, "tⁿ"), (5, "tⁿ+½"), (10, "tⁿ⁺¹")]:
           a.plot([x, x], [-0.95, -0.8], color=C["grey"], lw=1)
           text(a, x, -1.2, lab, color=C["grey"])
       a.plot([0, 10], [-0.875, -0.875], color=C["light"], lw=1)
       mark(a, "primal", 0, 1, color=C["grey"]); text(a, 0, 0.65, "Bⁿ", color=C["grey"])
       arrow(a, (0.3, 1), (9.6, 1), color=C["B"], lw=2.2)
       mark(a, "primal", 10, 1, color=C["B"]); text(a, 10, 1.45, flab, color=C["B"], bold=True, ha="right")
       text(a, 0.3, 1.25, far, color=C["B"], ha="left", va="bottom", fs=10)
       if push:
           mark(a, "dual", 5, 1, color=C["E"]); text(a, 5, 0.6, f"average {k + 1} of tⁿ and tⁿ⁺¹ fields", color=C["E"])
           mark(a, "primal", 0, 0, color=C["ion"])
           arrow(a, (0.3, 0), (9.6, 0), color=C["ion"], lw=2.2)
           mark(a, "primal", 10, 0, color=C["ion"])
           text(a, 10, -0.4, push, color=C["ion"], ha="right")
           text(a, 0.3, -0.4, "push with the average", color=C["ion"], ha="left")
       a.set_xlim(-0.6, 10.6); a.set_ylim(-1.45, 1.8)
   finish(fig, ax)

In equations, with :math:`\mathcal{O}(\mathbf{B}, \text{moments})` the electric field from
Ohm's law (:doc:`/theory/hybridpic`):

.. math::

   \begin{aligned}
   \mathbf{B}^{*} &= \mathbf{B}^n - \Delta t\,\nabla\times\mathbf{E}^n, &
   \mathbf{E}^{*} &= \mathcal{O}(\mathbf{B}^{*}, n^n, \mathbf{F}^n) \\
   \text{push with } &\tfrac12(\mathbf{E}^n + \mathbf{E}^{*}),\ \tfrac12(\mathbf{B}^n + \mathbf{B}^{*})
   & &\rightarrow\ \tilde n^{n+1}, \tilde{\mathbf{F}}^{n+1} \\
   \mathbf{B}^{**} &= \mathbf{B}^n - \tfrac{\Delta t}{2}\nabla\times(\mathbf{E}^n + \mathbf{E}^{*}), &
   \mathbf{E}^{**} &= \mathcal{O}(\mathbf{B}^{**}, \tilde n^{n+1}, \tilde{\mathbf{F}}^{n+1}) \\
   \text{push with } &\tfrac12(\mathbf{E}^n + \mathbf{E}^{**}),\ \tfrac12(\mathbf{B}^n + \mathbf{B}^{**})
   & &\rightarrow\ n^{n+1}, \mathbf{F}^{n+1} \\
   \mathbf{B}^{n+1} &= \mathbf{B}^n - \tfrac{\Delta t}{2}\nabla\times(\mathbf{E}^n + \mathbf{E}^{**}), &
   \mathbf{E}^{n+1} &= \mathcal{O}(\mathbf{B}^{n+1}, n^{n+1}, \mathbf{F}^{n+1})
   \end{aligned}

Here :math:`n` is the ion charge density and :math:`\mathbf{F}` the ion charge flux. The
first push only serves to predict the moments; its particles are thrown away. The second
push is the real one.

The MHD integrators
-------------------

The MHD solver computes, from a state :math:`U` (density, momentum, total energy, B), the
change rate :math:`L(U)` given by the finite-volume fluxes and the electric field on the
edges (:doc:`mhd_finite_volume`). One *Euler stage* is
:math:`U + \Delta t\,L(U)`. The integrators chain Euler stages. You choose one with
``mhd_timestepper`` (required for MHD). In the names, RK means Runge–Kutta; TVD (total-variation
diminishing) and SSP (strong-stability preserving) mean that the scheme adds no new
oscillations when a single Euler stage does not.

.. list-table::
   :header-rows: 1

   * - ``mhd_timestepper``
     - Stages
     - Scheme
   * - ``"Euler"``
     - 1
     - :math:`U^{n+1} = U^n + \Delta t\,L(U^n)` (first order)
   * - ``"TVDRK2"``
     - 2
     - :math:`U^{(1)} = U^n + \Delta t\,L(U^n)`;
       :math:`U^{n+1} = U^n + \tfrac{\Delta t}{2}\big(L(U^n) + L(U^{(1)})\big)` (second-order Runge–Kutta)
   * - ``"TVDRK3"``
     - 3
     - :math:`U^{(1)} = U^n + \Delta t\,L(U^n)`;
       :math:`U^{(2)} = \tfrac34 U^n + \tfrac14\big(U^{(1)} + \Delta t\,L(U^{(1)})\big)`;
       :math:`U^{n+1} = \tfrac13 U^n + \tfrac23\big(U^{(2)} + \Delta t\,L(U^{(2)})\big)`
       (third-order strong-stability-preserving Runge–Kutta)
   * - ``"SSPRK4_5"``
     - 5
     - five-stage, fourth-order strong-stability-preserving Runge–Kutta;
       coefficients in
       :src:`src/amr/solvers/time_integrator/ssprk4_5_integrator.hpp#L129-L153`

The name is not case-sensitive. In every scheme, the magnetic field is advanced with the
same weighted sum of edge electric fields as the other variables use for their fluxes, so
its divergence is not changed.

Choosing a stable time step
---------------------------

The limit is a :term:`CFL condition`.

**Hybrid.** Only a fixed time step is available: asking for ``"adaptive"`` with a hybrid
model stops the run during initialization. The limit comes
from the fastest wave, the whistler, whose frequency grows like :math:`k^2`. At the grid
scale this gives :math:`\Delta t \propto \Delta x^2`: halving the cell size requires a time
step four times smaller.

**MHD, adaptive.** At every step (and once at start-up), PHARE computes on each level the
smallest, over the level's cells, of

.. math::

   \Delta t_{\text{wave}} = \frac{C_{\text{wave}}}
   {\sum_d \big(|v_d| + c_{f,d} + c_{w,d}\big)/\Delta x_d},
   \qquad
   \Delta t_{\text{diff}} = \frac{C_{\text{diff}}}{2\eta\sum_d 1/\Delta x_d^2},

where :math:`c_{f,d}` is the fast magnetosonic speed along direction :math:`d`,
:math:`c_{w,d}` a whistler speed at the grid scale (only with ``hall=True``), and the
diffusive limit only applies when the resistivity :math:`\eta` is not zero
(:src:`src/amr/solvers/solver_mhd.hpp#L386-L459`).
:math:`C_{\text{wave}}` and :math:`C_{\text{diff}}` are ``cfl_wave`` and ``cfl_diffusive``,
between 0 and 1. Each level's limit is then converted to a coarsest-level time step by
multiplying it by :math:`4^{\ell}` (level :math:`\ell` takes :math:`4^{\ell}` steps per
coarsest step, see below), and the smallest over all levels and all MPI ranks is used
(:src:`src/amr/multiphysics_integrator.hpp#L435-L461`).

Time steps on refined levels
----------------------------

A refined level has cells half as large (refinement ratio 2). PHARE divides the time step
by the *square* of the ratio, so each finer level takes four steps for one step of the
level above (:src:`src/amr/multiphysics_integrator.hpp#L420-L428`).
The reason is the whistler limit above, :math:`\Delta t \propto \Delta x^2`. How the levels
are stepped and kept in sync is on :doc:`amr`.

Options
-------

.. list-table::
   :header-rows: 1

   * - Option
     - Meaning
     - Default
   * - ``time_step``
     - time step of the coarsest level, or a dict ``{"mode": "constant" | "adaptive", ...}``
     - none
   * - ``time_step_nbr``, ``final_time``
     - number of steps, end time (give two of the three)
     - none
   * - ``"cfl_wave"`` (key of the ``time_step`` dict)
     - adaptive mode: safety factor on the wave limit, > 0
     - none (required in adaptive mode)
   * - ``"cfl_diffusive"`` (key of the ``time_step`` dict)
     - adaptive mode: safety factor on the diffusion limit
     - ``cfl_wave``
   * - ``mhd_timestepper``
     - MHD integrator: ``Euler``, ``TVDRK2``, ``TVDRK3``, ``SSPRK4_5``
     - none (required for MHD)

Common mistakes
---------------

* **Adaptive time step in a hybrid run** — the run stops during initialization with
  ``Adaptive time-stepping not ready in hybrid.`` — use a fixed ``time_step``.
* **Adaptive time step with** ``time_step_nbr`` — rejected; give ``final_time`` instead.
* **Adaptive time step and diagnostics** — there is no fixed output grid; give explicit
  ``write_timestamps`` to each diagnostic (:doc:`/user_guide/diagnostics`).

.. seealso:: How this works inside: :doc:`/dev/one_time_step`.
