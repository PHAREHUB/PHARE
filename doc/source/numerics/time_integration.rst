=========
Time step
=========

**Goal:** set the time step of a hybrid or MHD run and keep it stable. How each model
advances over one step is on :doc:`hybrid_time_integration` and :doc:`mhd_time_integration`.

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

Hybrid and MHD integrators
--------------------------

* **Hybrid:** a predictor-predictor-corrector scheme (:term:`PPC`), on
  :doc:`hybrid_time_integration`.
* **MHD:** a Runge–Kutta integrator chosen with ``mhd_timestepper`` (required for MHD):
  ``Euler``, ``TVDRK2``, ``TVDRK3`` or ``SSPRK4_5``, on :doc:`mhd_time_integration`.

Choosing a stable time step
---------------------------

The limit is a :term:`CFL condition`.

**Hybrid.** Only a fixed time step is available: asking for ``"adaptive"`` with a hybrid
model stops the run during initialization. The limit comes
from the fastest wave, the whistler, whose frequency grows like :math:`k^2`. At the grid
scale this gives :math:`\Delta t \propto \Delta x^2`: halving the cell size requires a time
step four times smaller.

**MHD, adaptive.** At every step PHARE computes, on each level, the smallest over the
cells of a wave limit (fast magnetosonic speed, flow speed and, with Hall, a whistler speed
at the grid scale) and, when the resistivity is not zero, a diffusive limit. The formulas
are on :doc:`mhd_time_integration`. Each level's limit is then converted to a
coarsest-level time step by multiplying it by :math:`4^{\ell}` (level :math:`\ell` takes
:math:`4^{\ell}` steps per coarsest step, see below), and the smallest over all levels and
all MPI ranks is used (:src:`src/amr/multiphysics_integrator.hpp#L477-L503`).

Time steps on refined levels
----------------------------

A refined level has cells half as large (refinement ratio 2). PHARE divides the time step
by the *square* of the ratio, so each finer level takes four steps for one step of the
level above (:src:`src/amr/multiphysics_integrator.hpp#L462-L470`).
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
