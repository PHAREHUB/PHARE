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

**MHD, adaptive.** The stable step of a level is the smallest over its cells and over
all MPI ranks of a wave limit (fast magnetosonic speed, flow speed and, with Hall, a
whistler speed at the grid scale) and, when the resistivity is not zero, a diffusive limit.
The formulas are on :doc:`mhd_time_integration`. The coarsest level computes its stable
step before each of its steps and takes it as its step. Each finer level computes its own
after each of its sub-steps and uses it to size the next ones (below, and
:doc:`/dev/subcycling`) (:pr:`1570`).

Time steps on refined levels
----------------------------

A refined level has cells half as large (refinement ratio 2). With a fixed time step PHARE
divides the time step by the *square* of the ratio, so each finer level takes four steps for
one step of the level above (:src:`src/amr/multiphysics_integrator.hpp#L420-L428`).
The reason is the whistler limit above, :math:`\Delta t \propto \Delta x^2`. With an
adaptive time step each finer level splits the step of the level above into equal
sub-steps no longer than its own stable step: about two in ideal MHD (limit
:math:`\propto \Delta x`), up to four with Hall. How the levels are stepped and kept in sync
is on :doc:`amr`.

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
   * - ``"fine_level_substeps_max"`` (key of the ``time_step`` dict)
     - adaptive mode: stop the run when the stable step of a level is more than this many
       times smaller than that of the level above; integer ≥ 1
     - 1000
   * - ``"level0_dt_shrink_max"`` (key of the ``time_step`` dict)
     - adaptive mode: stop the run when the stable step of the coarsest level is more than
       this many times smaller than its first one (at the start or restart time); integer ≥ 1
     - 1000
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
* **Adaptive time step that collapses** — the run stops with ``stable dt on level ...``
  naming ``fine_level_substeps_max`` or ``level0_dt_shrink_max``, or with
  ``non-finite wave speed on level ...`` (the log names the first such cell). These stop a
  run whose step is shrinking toward zero; raise the limits only if such a small step is
  expected.

.. seealso:: How this works inside: :doc:`/dev/one_time_step`.
