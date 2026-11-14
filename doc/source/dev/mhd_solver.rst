.. _dev-mhd-solver:

====================
MHD solver internals
====================

.. container:: breadcrumb

   You are here: :doc:`/dev/index` › :doc:`/dev/one_time_step` › MHD solver internals

**In one sentence:** ``SolverMHD`` advances one :term:`level` of the MHD state with a Godunov
finite-volume scheme for density, momentum and energy, and constrained transport for
the magnetic field, wrapped in a Runge-Kutta time integrator chosen at run time.

Terms: :term:`Godunov method`, :term:`Riemann solver`, :term:`constrained transport`,
:term:`reflux`.

The numerical method itself (reconstruction, Riemann solvers, constrained transport,
the formulas of each integrator) is explained for users on
:doc:`/numerics/mhd_finite_volume` and :doc:`/numerics/time_integration`. This page is
about how the code is put together.

Why it is built this way
------------------------

A Godunov MHD step has many interchangeable pieces: the reconstruction, its slope
limiter, the Riemann solver, the Hall term, dissipation, the time integrator. Some are
inside the innermost loop over faces; making them runtime choices would cost a branch
per face. So the inner pieces are template parameters, fixed at compile time, and only
the outer choices (time integrator, resistivity values) are read at run time.

The picture
-----------

One call to ``SolverMHD::advanceLevel``, for one level:

.. mermaid::

   flowchart TD
     A["SolverMHD::advanceLevel"] --> T["TimeIntegrator<br/>Euler, TVDRK2, TVDRK3, SSPRK4_5"]
     T -->|"each stage"| E["Euler stage"]
     E --> C["ComputeFluxes"]
     C --> C1["to primitive<br/>Ampère (if needed)<br/>dissipative E (if needed)"]
     C1 --> G["Godunov fluxes<br/>reconstruct, Riemann,<br/>save CT coefficients"]
     G --> CT["back to conservative<br/>constrained transport: E on edges<br/>fill E ghosts"]
     CT --> U["EulerUsingComputedFlux<br/>finite-volume update, Faraday,<br/>fill state ghosts"]
     click E href "mhd_solver.html#one-euler-stage" "One Euler stage"
     click CT href "messengers.html" "Messengers"
     click U href "messengers.html" "Messengers"

How it works
------------

#. **Before the step.** ``SolverMHD::prepareStep`` copies ``rho``, ``rhoV``, ``B`` and
   ``Etot`` into ``stateOld_`` and records the time
   (:src:`src/amr/solvers/solver_mhd.hpp#L276-L306`). This copy is only used later to
   *reflux* (see below).
#. **The step.** ``SolverMHD::advanceLevel`` calls the time integrator once
   (:src:`src/amr/solvers/solver_mhd.hpp#L360-L383`). The integrator is chosen at run time
   from ``mhd_timestepper`` by a small factory
   (:src:`src/amr/solvers/time_integrator/time_integrator.hpp#L25-L47`).
#. **Each integrator stage is one Euler step** (``Euler::operator()``,
   :src:`src/amr/solvers/time_integrator/euler.hpp#L23-L32`): compute the fluxes and the
   edge electric field from a state, then apply them over a time step. Multi-stage
   integrators chain Euler stages and combine states.
#. **After the step.** The flux and electric field that produced the step are summed
   over the fine sub-steps and later used to correct the coarser level (refluxing,
   :doc:`/dev/coarsening`).

One Euler stage
^^^^^^^^^^^^^^^

``ComputeFluxes::operator()``
(:src:`src/amr/solvers/time_integrator/compute_fluxes.hpp#L63-L90`), in this order:

#. Convert the conservative variables (``rho``, ``rhoV``, ``Etot``) to primitive ones
   (``V``, ``P``) on the whole ghost box.
#. If the Hall term, resistivity or hyper-resistivity is on, compute
   :math:`\mathbf{J} = \nabla\times\mathbf{B}` (Ampère) and, if dissipative, the
   dissipative electric field :math:`\eta\mathbf{J} - \nu\nabla^2\mathbf{J}` on the edges.
#. ``Godunov``: for each direction, on each face, reconstruct the left and right states,
   compute the physical fluxes, call the Riemann solver, and save what constrained
   transport needs (:src:`src/core/numerics/godunov_fluxes/godunov_fluxes.hpp#L99-L159`).
#. Convert back to conservative variables.
#. Constrained transport: build E on the cell edges
   (``UpwindConstrainedTransport``,
   :src:`src/core/numerics/constrained_transport/upwind_constrained_transport.hpp#L34-L53`),
   then fill E's ghosts through the messenger.

Then ``EulerUsingComputedFlux::operator()``
(:src:`src/amr/solvers/time_integrator/euler_using_computed_flux.hpp#L24-L35`):

#. Finite-volume update of ``rho``, ``rhoV``, ``Etot``:
   :math:`U^{new}_i = U_i - \frac{\Delta t}{\Delta x}(F_{i+1/2} - F_{i-1/2})` per direction.
#. Faraday on the faces: :math:`\mathbf{B}^{new} = \mathbf{B} - \Delta t\,\nabla\times\mathbf{E}`.
#. Fill all ghosts of the new state through the messenger.

The tricks
----------

Compile-time scheme, run-time integrator
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** Reconstruction, limiter, Riemann solver and Hall term sit in the loop over
faces. Choosing them at run time costs a branch per face per stage.

**Fix.** They are template parameters, resolved from enums by ``mhd_resolver.hpp``
(:src:`src/amr/solvers/mhd_resolver.hpp#L35-L131`). Python picks the matching
pre-compiled module, named after the combination
(:src:`pyphare/pyphare/cpp/__init__.py#L8-L35`); only the combinations listed in
:src:`res/sim/all.txt` are built by default. The time integrator, resistivity and
hyper-resistivity are runtime choices: the integrator is a virtual class behind a
factory, and the ``eta > 0`` / ``nu > 0`` flags are turned into compile-time branches
once per call (``Constexprifier``, as in Ohm's law).

**History.** :pr:`1233` (runtime integrator and resistivity).

Invalid reconstruction/limiter pairs do not compile
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** Only ``Linear`` reconstruction uses a slope limiter. A half-configured
build (say, a reconstruction with no limiter set) should not silently pick something.

**Fix.** ``SlopeLimiterSelector`` is only declared; every valid pair is an explicit
specialisation (``Linear`` with ``VanLeer`` or ``MinMod``; the others with ``None``).
Any other pair fails to compile (:src:`src/amr/solvers/mhd_resolver.hpp#L70-L110`).

The normal B is not reconstructed
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** The Riemann problem at an x-face needs left and right values of all
variables, including :math:`B_x`. But :math:`B_x` already lives on that face, and it
must be the same on both sides for the MHD Riemann problem to make sense.

**Fix.** The face value of the normal component is used as is for both sides; only the
transverse components are averaged to cell centres and then reconstructed
(:src:`src/core/numerics/reconstructions/reconstructor.hpp#L64-L95`).

E on the edges comes from the Riemann solver's coefficients
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** Faraday needs E on the cell edges. Averaging the face fluxes of B to the
edges is simple but not upwind, and adds oscillations.

**Fix.** Upwind constrained transport (UCT). On each face, the Riemann solver exposes four
weights :math:`(a_L, a_R, d_L, d_R)` and the averaged transverse velocity; ``Godunov``
saves them per face (:src:`src/core/numerics/godunov_fluxes/godunov_fluxes.hpp#L129-L146`).
``UpwindConstrainedTransport`` then assembles each edge E from the saved values of the
faces that meet at that edge, with B and the velocity reconstructed again at the edge.
The HLL weights follow the Idefix code (a comment says so,
:src:`src/core/numerics/riemann_solvers/hll.hpp#L173-L191`). Because B changes only by the
discrete curl of an edge E, its discrete divergence is not changed by the update.

The stage fluxes are kept for refluxing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** After the fine level has advanced, the coarse level must be corrected with
the fine fluxes through the coarse-fine border. With a multi-stage integrator, "the flux of
the step" is a weighted sum of stage fluxes.

**Fix.** Every multi-stage integrator accumulates its stage fluxes and edge E with the Runge-Kutta
weights (``accumulateButcherFluxes_``,
:src:`src/amr/solvers/time_integrator/base_mhd_timestepper.hpp#L113-L122`) and finishes with
a single Euler update from :math:`U^n` using that sum. The same sum is what the level
hands over for refluxing, so the coarse correction is consistent with the step taken.
``SolverMHD::reflux`` redoes that Euler update from ``stateOld_`` with the corrected fluxes
(:src:`src/amr/solvers/solver_mhd.hpp#L346-L357`).

Intermediate stages get their ghosts at their own time
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** Combining two states (for example :math:`U^{(2)} = \tfrac34 U^n + \tfrac14 U^{(1)}`
in TVDRK3) only touches the domain values, and that stage lives at :math:`t + \Delta t/2`.

**Fix.** The integrator stamps the combined state with its own time and fills its ghosts
there, so that coarse-fine ghosts are interpolated in time at the stage's abscissa
(:src:`src/amr/solvers/time_integrator/tvdrk3_integrator.hpp#L51-L59`).

Ghost width follows the reconstruction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** The reconstruction reads several cells on each side of a face, and it also
reads J, which Ampère computes on the ghost box minus one layer.

**Fix.** The ghost width is the stencil width plus one, rounded up to an even number
(:src:`src/core/utilities/ghost_width_calculator.hpp#L47-L68`), and ``Godunov`` checks it at
compile time (:src:`src/core/numerics/godunov_fluxes/godunov_fluxes.hpp#L85-L86`). See
:doc:`/dev/gridlayout`.

What the solver does not do
^^^^^^^^^^^^^^^^^^^^^^^^^^^

* No density or pressure floor: a pressure clamp is present but commented out
  (:src:`src/core/numerics/primite_conservative_converter/to_primitive_converter.hpp#L24-L32`).
* The wave speeds in the Riemann solvers do not include the whistler speed (it is set to
  zero, :src:`src/core/numerics/riemann_solvers/rusanov.hpp#L123-L130`); the whistler speed
  only enters the time-step estimate when the Hall term is on.
* Hyper-resistivity requires the Hall term: ``nu > 0`` without ``hall=True`` throws
  (:src:`src/amr/solvers/time_integrator/compute_fluxes.hpp#L56-L60`).

Go deeper
---------

* :doc:`/dev/gridlayout`: the MHD centering table (cell-centred fluid, B on faces, E on
  edges).
* :doc:`/numerics/mhd_finite_volume`: the scheme and its options.

Related
-------

* :doc:`/dev/one_time_step` (up): where the solver is called.
* :doc:`/dev/messengers`: how MHD ghosts are filled (``fillElectricGhosts``,
  ``fillStateGhosts``).
* :doc:`/dev/field_refinement`: refining B so that its divergence stays zero.
* :doc:`/dev/coarsening`: refluxing.
* :doc:`/dev/electromagnetic`: the hybrid counterpart.
