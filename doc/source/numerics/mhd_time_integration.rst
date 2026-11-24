==========================================
MHD time integration: stages and time step
==========================================

**Goal:** know how the MHD solver advances one time step, what each of the four
integrators costs and gives, and how the time step is chosen.

**In one sentence:** the spatial scheme (:doc:`mhd_finite_volume`) gives the rate of
change of the state; a Runge–Kutta integrator chains several Euler stages built on it, to
make the steps more accurate in time without making them less stable.

.. contents:: On this page
   :local:
   :depth: 1

One Euler stage
---------------

Call :math:`L(U)` the rate of change of the state :math:`U` (density, momentum, total
energy, magnetic field) given by the fluxes through the faces and the electric field on the
edges. One *Euler stage* computes the fluxes from the current state, then updates

* the fluid variables with the divergence of the fluxes,
* B with the curl of the edge electric field, :math:`\mathbf B \leftarrow \mathbf B -
  \Delta t\,\nabla\times\mathbf E`, where E comes from the same stage

(:src:`src/amr/solvers/time_integrator/euler_using_computed_flux.hpp#L24-L35`). All the
integrators below chain such stages. The cost of a step is the number of stages, because
each needs one full flux computation.

The four integrators
--------------------

``mhd_timestepper`` takes ``"Euler"``, ``"TVDRK2"``, ``"TVDRK3"`` or ``"SSPRK4_5"``. The name
is not case sensitive; an unknown name stops the run with an error
(:src:`pyphare/pyphare/pharein/initialize/dict_populator.py#L50-L61`). "RK" stands for
Runge–Kutta; "TVD" and "SSP" (strong-stability-preserving) mean that, when a single Euler
stage adds no new oscillation, the whole step adds none either, provided the time step is
small enough.

.. list-table::
   :header-rows: 1
   :widths: 40 30 30

   * - ``mhd_timestepper``
     - Flux computations per step
     - Order in time
   * - ``Euler``
     - 1
     - 1
   * - ``TVDRK2``
     - 2
     - 2
   * - ``TVDRK3``
     - 3
     - 3
   * - ``SSPRK4_5``
     - 5
     - 4

With :math:`U^n` the state at the start of the step:

* ``Euler``: :math:`U^{n+1}=U^n+\Delta t\,L(U^n)`
  (:src:`src/amr/solvers/time_integrator/euler_integrator.hpp#L20-L35`).
* ``TVDRK2``: :math:`U^{(1)}=U^n+\Delta t\,L(U^n)`, then
  :math:`U^{n+1}=U^n+\tfrac{\Delta t}{2}\big(L(U^n)+L(U^{(1)})\big)`
  (:src:`src/amr/solvers/time_integrator/tvdrk2_integrator.hpp#L24-L44`).
* ``TVDRK3``: :math:`U^{(1)}=U^n+\Delta t\,L(U^n)`;
  :math:`U^{(2)}=\tfrac34U^n+\tfrac14\big(U^{(1)}+\Delta t\,L(U^{(1)})\big)`;
  :math:`U^{n+1}=\tfrac13U^n+\tfrac23\big(U^{(2)}+\Delta t\,L(U^{(2)})\big)`
  (:src:`src/amr/solvers/time_integrator/tvdrk3_integrator.hpp#L38-L76`).
* ``SSPRK4_5``: five stages of the form
  :math:`U^{(k)}=a\,U^n+b\,U^{(k-1)}+c\,\Delta t\,L(U^{(k-1)})`, with the fifteen
  coefficients listed in the source
  (:src:`src/amr/solvers/time_integrator/ssprk4_5_integrator.hpp#L130-L144`). The code
  comment names it the five-stage, fourth-order scheme of Spiteri and Ruuth
  (:src:`src/amr/solvers/time_integrator/ssprk4_5_integrator.hpp#L146-L149`).

In every case the step ends as :math:`U^{n+1}=U^n+\Delta t\sum_k \beta_k\,L(U^{(k)})`: the
rates of change of all stages are added with fixed weights :math:`\beta_k`. The figure
shows these weights, each at the time :math:`t+c_k\Delta t` that its stage represents.

.. plot::
   :include-source: False
   :alt: For each integrator, the weight of the rate of change of each stage in the final update, plotted at the time of the stage.

   import numpy as np
   import matplotlib.pyplot as plt
   from phare_draw import C

   # SSPRK4_5 coefficients: ssprk4_5_integrator.hpp:130-144; abscissae :150-153
   w0=0.391752226571890; w11=0.555629506348765; w12=0.368410593050371
   w21=0.379898148511597; w22=0.251891774271694; w31=0.821920045606868; w32=0.544974750228521
   c1=w0; c2=w11*c1+w12; c3=w21*c2+w22; c4=w31*c3+w32
   beta5=[0.1468118760847865, 0.24848290944497606, 0.10425883033198079, 0.27443890090135015, 0.226007483236906]
   rows = [("Euler", [0.0], [1.0]),
           ("TVDRK2", [0.0, 1.0], [0.5, 0.5]),
           ("TVDRK3", [0.0, 1.0, 0.5], [1/6, 1/6, 2/3]),
           ("SSPRK4_5", [0.0, c1, c2, c3, c4], beta5)]
   fig, axes = plt.subplots(4, 1, figsize=(3.4, 6.4), sharex=True, sharey=True)
   for ax, (name, c, b) in zip(axes, rows):
       ax.axhline(0, color=C["grey"], lw=1.2)
       for ck, bk in zip(c, b):
           ax.plot([ck, ck], [0, bk], color=C["L1"], lw=3)
           ax.plot([ck], [bk], "o", color=C["L1"], ms=5)
           off = 0.04 + (0.17 if name == "SSPRK4_5" and ck == c3 else 0.0)
           ax.text(ck, bk + off, f"{bk:.3f}".rstrip("0").rstrip("."), ha="center", fontsize=7)
       ax.set_title(name, fontsize=10, loc="left")
       ax.set_ylim(0, 1.25); ax.set_yticks([])
       for s in ("top", "right", "left"):
           ax.spines[s].set_visible(False)
       ax.set_xticks([0, 0.5, 1]); ax.set_xticklabels(["t", "t + Δt/2", "t + Δt"], fontsize=8)
   fig.tight_layout()

* **Weights.** The weights add up to 1. For ``TVDRK3`` they are 1/6, 1/6, 2/3 for
  :math:`U^n`, :math:`U^{(1)}`, :math:`U^{(2)}` (:src:`src/amr/solvers/time_integrator/tvdrk3_integrator.hpp#L94-L97`).
  For ``SSPRK4_5`` they are, in the same order, 0.1468, 0.2485, 0.1043, 0.2744 and 0.2260.
* **Times.** Each stage value stands for a time :math:`t+c_k\Delta t`, the time at which it
  fills its ghost cells on refined levels. For ``TVDRK3`` the second stage is at
  :math:`t+\Delta t/2`; for ``SSPRK4_5`` the four stages are at :math:`c_k` = 0.392, 0.586,
  0.475 and 0.935 (:src:`src/amr/solvers/time_integrator/tvdrk3_integrator.hpp#L64`,
  :src:`src/amr/solvers/time_integrator/ssprk4_5_integrator.hpp#L150-L153`).

Order in time
-------------

The order of an integrator is how fast its error falls when the time step is halved:
:math:`2^p` for order :math:`p`. The orders in the table follow from the coefficients: for
``SSPRK4_5``, the script ``doc/source/validation/scripts/ssprk_check.py`` rebuilds the
five stages from the coded coefficients and finds the expected fourth order on a nonlinear
test problem (4.00).

A *documentation run*, not part of the test suite, measures them on the solver itself. It
uses the small Alfvén wave of :doc:`mhd_reconstruction` on 64 cells over ten periods, with
the ``MP5`` reconstruction, and halves the time step:

* ``Euler``: the error falls by 2.2, 2.1 and 2.0 for each halving of :math:`\Delta t`
  (order 1);
* ``TVDRK2``: by 4.0 each time (order 2);
* ``TVDRK3``: by 8 each time once the space error is removed (order 3; the space error is taken from the ``SSPRK4_5`` run with the same step, as in ``validation/scripts/mhd_time_order_ratios.py``);
* ``SSPRK4_5``: the time error is below the space error, so its order cannot be read there.

In this run ``Euler`` gave NaN for time steps from 0.003125 and ``TVDRK2`` at 0.0125. The
fast speed is 1 there, so these are Courant numbers of 0.2 and 0.8. This is what was
observed in this run, not a stability limit.

The error of a run is dominated by its lowest order, in space or in time. That is why the
convergence test of :doc:`mhd_reconstruction` uses ``SSPRK4_5`` with a small time step: the
time error stays below the space error.

Divergence of B
---------------

B is advanced with the electric field of each stage, added with the same weights
:math:`\beta_k` as the fluxes (:src:`src/amr/solvers/time_integrator/base_mhd_timestepper.hpp#L115-L124`),
so the integrator does not change :math:`\nabla\cdot\mathbf B`. Why, and what happens
across refined levels: :doc:`div_b`.

Memory
------

An integrator keeps copies of the state while it works: none for ``Euler``, 1 for
``TVDRK2``, 2 for ``TVDRK3``, 4 for ``SSPRK4_5``
(:src:`src/amr/solvers/time_integrator/euler_integrator.hpp#L21`,
:src:`src/amr/solvers/time_integrator/tvdrk2_integrator.hpp#L19`,
:src:`src/amr/solvers/time_integrator/tvdrk3_integrator.hpp#L30`,
:src:`src/amr/solvers/time_integrator/ssprk4_5_integrator.hpp#L30`).

Choosing the time step
----------------------

A fixed time step is given as in the hybrid case (:doc:`time_integration`). The MHD solver
also offers an *adaptive* step, recomputed at every step (and once at start-up)::

   time_step={"mode": "adaptive", "cfl_wave": 0.8}

With refined levels, each finer level sizes its sub-steps from its own stable step, and the
coarse level is corrected with the fine fluxes weighted by each sub-step's share of the
coarse step, so mass, momentum and energy stay conserved across levels (:doc:`amr`).

On each level, PHARE computes the smallest over the cells of

.. math::

   \Delta t_{\text{wave}} = \frac{C_{\text{wave}}}
   {\sum_d \big(|v_d| + c_{f,d} + c_{w,d}\big)/\Delta x_d},
   \qquad
   \Delta t_{\text{diff}} = \frac{C_{\text{diff}}}{2\eta\sum_d 1/\Delta x_d^2},

where :math:`c_{f,d}` is the fast magnetosonic speed along direction :math:`d` and
:math:`c_{w,d}` what the whistler adds to the Alfvén speed at the grid scale, only with
``hall=True``:

.. math::

   c_{w,d} = v_A\Big(\sqrt{1+x_d^2}+x_d-1\Big),\qquad
   x_d = \frac{\pi}{2\,\Delta x_d\sqrt{\rho}},\qquad
   v_A = \frac{|\mathbf B|}{\sqrt{\rho}}.

:math:`v_A(\sqrt{1+x^2}+x)` is the phase speed of a whistler of wavenumber
:math:`\pi/\Delta x`, the ion inertial length being :math:`1/\sqrt{\rho}`; the Alfvén part
is already in :math:`c_{f,d}`. On fine grids :math:`x_d` is large, :math:`c_{w,d}` grows
like :math:`1/\Delta x` and the step goes as :math:`\Delta x^2` (:pr:`1570`). The
diffusive limit only applies when ``eta`` is not zero; hyper-resistivity gives no limit
(:src:`src/amr/solvers/solver_mhd.hpp#L386-L459`,
:src:`src/core/numerics/riemann_solvers/mhd_speeds.hpp#L23-L27`).
:math:`C_{\text{wave}}` and :math:`C_{\text{diff}}` are ``cfl_wave`` and ``cfl_diffusive``,
positive numbers (the docstring advises values in (0, 1], 1 being the estimated limit; ``cfl_diffusive`` defaults to
``cfl_wave``) (:src:`pyphare/pyphare/pharein/simulation.py#L889-L900`). The estimate is not a
guarantee: in the documentation run above the first- and second-order integrators failed at Courant
numbers of 0.2 (``Euler``) and 0.8 (``TVDRK2``). How the coarsest level's limit becomes
its step and the finer levels' limits set their sub-steps is on
:doc:`time_integration` and :doc:`/dev/subcycling`; how the fine-level flux sums are used to
correct the coarse level is on :doc:`/dev/refluxing`.

How to choose
-------------

* **Cost.** More stages cost more flux computations per step and more memory (table and
  section above). The scheme is second order at best in space (:doc:`mhd_finite_volume`),
  so an integrator of order 3 or 4 does not raise the order of the scheme beyond that.
* **With** ``WENOZ`` **or** ``MP5``: the literature recommends an integrator of order three
  or more for WENO-Z (Wang and Spiteri 2007). PHARE does not enforce it. The tests that
  use ``WENOZ`` use ``TVDRK3`` or ``SSPRK4_5``; the shock-tube and Orszag–Tang tests use
  ``TVDRK2`` with ``Linear``.
* **Step size.** ``cfl_wave`` and ``cfl_diffusive`` multiply the estimated limits above; the
  docstring advises values in (0, 1]. The first-order and second-order integrators failed at
  0.2 and 0.8 in the documentation run, so staying at or below 1 does not guarantee stability.
* **With** ``hall=True``: the discrete Hall operator linearised about a uniform B is skew (its
  spectrum is imaginary) and the upwind dissipation cannot damp the whistler at the grid scale once
  :math:`\Delta x \lesssim d_i`. ``Euler`` and ``TVDRK2``, which have no stability range on
  the imaginary axis, are then unstable whatever the step unless ``nu > 0``; use ``TVDRK3`` or
  ``SSPRK4_5``. PHARE does not check it. With the estimate above at :math:`k=\pi/\Delta x`,
  the linear stability limit of ``cfl_wave`` for the Hall term is :math:`\pi\rho_I/4`, with
  :math:`\rho_I` the reach of the integrator on the imaginary axis: 1.36 for ``TVDRK3``
  (:math:`\rho_I=\sqrt3`) and 2.57 for ``SSPRK4_5`` (:math:`\rho_I=3.28`).

.. seealso::
   :doc:`mhd_reconstruction`, :doc:`mhd_riemann_solvers`, :doc:`mhd_finite_volume`,
   :doc:`time_integration`.
