.. _dev-mhd-solver:

====================
MHD solver internals
====================

.. container:: breadcrumb

   You are here: :doc:`/dev/index` › :doc:`/dev/one_time_step` › MHD solver internals

**In one sentence:** ``SolverMHD`` advances one :term:`level` with a :term:`Godunov method`
for ρ, ρv, Etot and :term:`constrained transport` for B, inside a Runge–Kutta integrator
chosen at run time. The scheme for users: :doc:`/numerics/mhd_finite_volume`.

One Euler stage
---------------

.. plot::
   :include-source: False
   :caption: One stage on the MHD grid: fluid averages at cell centres (squares), B on
      faces (arrows). Look at what is computed on faces (panels 2–3), then on edges (4).

   from phare_draw import *
   from matplotlib.patches import Circle

   T = ["Cell averages ρ, ρv, Etot; Bx on faces",
        "Reconstruct left/right states at face i+½",
        "Riemann solver: flux F, UCT weights saved",
        "2D: Ez at a corner from its 4 faces",
        "Update: U −= Δt/Δx (F(i+½) − F(i−½)), B −= Δt curl E"]
   fig, ax = strip(5, titles=T)
   u = [1.0, 0.85, 0.45, 0.3]
   for k in (0, 1, 2, 4):
       a = ax[k]
       grid(a, (0, 0, 4, 1))
       for i, ui in enumerate(u):
           y = 0.25 + 0.5 * ui if k == 1 else 0.5   # panel 2: averages drawn at their value
           mark(a, "dual", i + 0.5, y, color=C["rho"])
       for x in (range(5) if k != 1 else [2]):
           arrow(a, (x - 0.18, 0.5), (x + 0.18, 0.5), color=C["B"], lw=2, ms=10)
       text(a, 0.5, 1.25, "U(i−1)", color=C["rho"]); text(a, 1.5, 1.25, "U(i)", color=C["rho"])
       text(a, 2.5, 1.25, "U(i+1)", color=C["rho"])
   a = ax[1]
   for i, (lo, hi) in enumerate([(1.05, 0.95), (0.95, 0.75), (0.6, 0.3)]):
       a.plot([i + 0.05, i + 0.95], [0.25 + 0.5 * lo, 0.25 + 0.5 * hi], color=C["rho"], lw=2)
   mark(a, "primal", 1.95, 0.25 + 0.5 * 0.75, color=C["rho"], size=7)
   mark(a, "primal", 2.05, 0.25 + 0.5 * 0.6, color=C["rho"], size=7)
   text(a, 1.7, 0.85, "U_L", color=C["rho"]); text(a, 2.3, 0.85, "U_R", color=C["rho"])
   text(a, 2, -0.3, "Bx(i+½) used as is on both sides", color=C["B"])
   a = ax[2]
   a.add_patch(Circle((2, 0.5), 0.16, fc=C["hi"], ec=C["ink"], lw=1.5, zorder=4))
   arrow(a, (1.6, 0.85), (2.4, 0.85), color=C["ink"], lw=2.4, ms=14)
   text(a, 2, 1.12, "F(i+½)", bold=True)
   text(a, 2, -0.3, "saved per face: aL, aR, dL, dR, transverse v", color=C["grey"])
   a = ax[3]
   grid(a, (0, 0, 2, 2))
   for y in (0.5, 1.5):
       arrow(a, (1 - 0.2, y), (1 + 0.2, y), color=C["B"], lw=2, ms=10)
   for x in (0.5, 1.5):
       arrow(a, (x, 1 - 0.2), (x, 1 + 0.2), color=C["B"], lw=2, ms=10)
   for (p, q) in [((1, 0.6), (1, 0.85)), ((1, 1.4), (1, 1.15)), ((0.6, 1), (0.85, 1)), ((1.4, 1), (1.15, 1))]:
       arrow(a, p, q, color=C["grey"], lw=1.4, ms=10)
   a.add_patch(Circle((1, 1), 0.12, fc=C["hi"], ec=C["E"], lw=2, zorder=4))
   text(a, 1.25, 1.25, "Ez", color=C["E"], bold=True)
   text(a, 1.3, 0.35, "x-faces", color=C["grey"]); text(a, 0.2, 1.3, "y-faces", color=C["grey"])
   a = ax[4]
   for x in (1, 2):
       arrow(a, (x - 0.3, 0.85), (x + 0.3, 0.85), color=C["ink"], lw=2, ms=12)
   a.add_patch(Circle((1.5, 0.5), 0.16, fc=C["hi"], ec=C["rho"], lw=2, zorder=4))
   text(a, 1, 1.12, "F(i−½)"); text(a, 2, 1.12, "F(i+½)")
   finish(fig, ax)

* ``ComputeFluxes`` does panels 1–4: to primitive on the ghost box; J (and the dissipative
  E) if Hall, η or ν; ``Godunov`` per direction; back to conservative; ``UpwindConstrainedTransport``;
  fill E ghosts (:src:`src/amr/solvers/time_integrator/compute_fluxes.hpp#L63-L90`).
* ``EulerUsingComputedFlux`` does panel 5 and fills the new state's ghosts
  (:src:`src/amr/solvers/time_integrator/euler_using_computed_flux.hpp#L24-L35`).
* Panel 2: the normal B is not reconstructed, only the transverse ones
  (:src:`src/core/numerics/reconstructions/reconstructor.hpp#L64-L95`).
* Panel 4: edge E from the saved weights; B (and, in 3D, the velocity) reconstructed again
  at the edge (:src:`src/core/numerics/constrained_transport/upwind_constrained_transport.hpp#L34-L53`).
  B changes only by a discrete curl, so its discrete divergence does not change.

The integrator around the stages
--------------------------------

.. plot::
   :include-source: False
   :caption: TVDRK3. Each stage's fluxes are added, with its weight, to one sum; the step is
      one Euler update from Uⁿ with that sum, and the same sum is kept for refluxing.

   from phare_draw import *

   fig, ax = strip(1, titles=["TVDRK3: three stages, one flux sum"], aspect=False, height=2.2)
   a = ax[0]
   a.plot([0, 10], [0, 0], color=C["grey"], lw=1.5)
   for x, lab in [(0, "tⁿ"), (5, "tⁿ + Δt/2"), (10, "tⁿ + Δt")]:
       a.plot([x, x], [-0.12, 0.12], color=C["grey"], lw=1.5)
       text(a, x, -0.4, lab)
   for x, lab, w in [(0, "F(Uⁿ)", "1/6"), (10, "F(U¹)", "1/6"), (5, "F(U²)", "2/3")]:
       mark(a, "primal", x, 0.6, color=C["rho"])
       text(a, x, 1.0, lab, color=C["rho"], bold=True)
       text(a, x, 1.4, "× " + w, color=C["ink"])
   text(a, 5, -0.85, "U² = ¾Uⁿ + ¼(U¹ + Δt F(U¹)): ghosts filled at tⁿ + Δt/2", color=C["grey"])
   a.set_xlim(-1.2, 11.2); a.set_ylim(-1.15, 1.7)
   finish(fig, ax)

* The integrator is chosen at run time from ``mhd_timestepper`` (``Euler``, ``TVDRK2``,
  ``TVDRK3``, ``SSPRK4_5``) by a factory
  (:src:`src/amr/solvers/time_integrator/time_integrator.hpp#L25-L47`).
* Each multi-stage integrator accumulates stage fluxes and E with its weights
  (``accumulateButcherFluxes_``, :src:`src/amr/solvers/time_integrator/base_mhd_timestepper.hpp#L113-L122`).
* ``SolverMHD::prepareStep`` copies the state to ``stateOld_``; ``reflux`` redoes the Euler
  update from it with the corrected fluxes (:src:`src/amr/solvers/solver_mhd.hpp#L276-L357`).

The tricks
----------

Compile-time scheme, run-time integrator
   **Problem.** Reconstruction, limiter, Riemann solver and Hall sit in the per-face loop.
   **Fix.** Template parameters (:src:`src/amr/solvers/mhd_resolver.hpp#L35-L131`); Python
   loads the compiled module for that combination (:src:`pyphare/pyphare/cpp/__init__.py#L8-L35`).
   Integrator, η, ν stay run-time. :pr:`1233` (runtime integrator and resistivity).

Invalid reconstruction/limiter pairs do not compile
   **Problem.** Only ``Linear`` uses a limiter; a half-configured build must not guess.
   **Fix.** Every valid pair is an explicit specialisation; any other fails to compile
   (:src:`src/amr/solvers/mhd_resolver.hpp#L70-L110`).

Stage ghosts at the stage's own time
   **Problem.** A combined state (U² in TVDRK3) has domain values only, and lives at t + Δt/2.
   **Fix.** It is stamped with its time and its ghosts filled there
   (:src:`src/amr/solvers/time_integrator/tvdrk3_integrator.hpp#L57-L66`).

Ghost width follows the reconstruction
   **Problem.** The stencil reads several cells, and J (computed on the ghost box minus one).
   **Fix.** Width = stencil + 1, rounded up to even, checked at compile time
   (:src:`src/core/utilities/ghost_width_calculator.hpp#L47-L68`,
   :src:`src/core/numerics/godunov_fluxes/godunov_fluxes.hpp#L85-L86`).

What the solver does not do
   No density or pressure floor (clamp commented out,
   :src:`src/core/numerics/primite_conservative_converter/to_primitive_converter.hpp#L24-L32`);
   no whistler speed in the Riemann solvers (:src:`src/core/numerics/riemann_solvers/rusanov.hpp#L123-L130`);
   ``nu > 0`` needs ``hall=True`` (:src:`src/amr/solvers/time_integrator/compute_fluxes.hpp#L56-L60`).

Go deeper
---------

* :doc:`/dev/gridlayout` — the MHD centering table.
* :doc:`/numerics/mhd_finite_volume` · :doc:`/numerics/time_integration`

Related
-------

* :doc:`/dev/one_time_step` (up) · :doc:`/dev/messengers` · :doc:`/dev/field_refinement` ·
  :doc:`/dev/coarsening` · :doc:`/dev/electromagnetic`
