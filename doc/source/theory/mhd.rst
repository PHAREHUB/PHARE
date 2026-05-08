
===================
The MHD formalism
===================

Magnetohydrodynamics (MHD) treats the plasma as a single conducting fluid. It is valid at scales
much larger than ion inertial lengths and gyroradii, where kinetic effects average out. In PHARE,
the MHD solver can run standalone or provide large-scale context on coarse AMR levels while the
Hybrid PIC solver handles kinetic physics on fine levels (see :doc:`amr`).


The ideal MHD equations
-----------------------

The ideal MHD equations express conservation of mass, momentum, energy, and magnetic flux. Written
in conservation form they read:

**Mass conservation**

.. math::

   \frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \mathbf{v}) = 0

**Momentum conservation**

.. math::

   \frac{\partial (\rho \mathbf{v})}{\partial t}
   + \nabla \cdot \left( \rho \mathbf{v}\mathbf{v} + P^* \mathbf{I}
     - \frac{\mathbf{B}\mathbf{B}}{\mu_0} \right) = 0

where :math:`P^* = P + \dfrac{B^2}{2\mu_0}` is the total (thermal + magnetic) pressure.

**Total energy conservation**

.. math::

   \frac{\partial E_\mathrm{tot}}{\partial t}
   + \nabla \cdot \left[ \left( E_\mathrm{tot} + P^* \right) \mathbf{v}
     - \frac{\mathbf{B}(\mathbf{v} \cdot \mathbf{B})}{\mu_0} \right] = 0

where the total energy density is

.. math::

   E_\mathrm{tot} = \frac{P}{\gamma - 1} + \frac{\rho v^2}{2} + \frac{B^2}{2\mu_0}

**Induction equation**

.. math::

   \frac{\partial \mathbf{B}}{\partial t} - \nabla \times (\mathbf{v} \times \mathbf{B}) = 0

Together these form a hyperbolic system of conservation laws that can be integrated with standard
finite-volume methods.


Finite volume discretization
-----------------------------

The MHD equations are integrated over each cell volume :math:`\Omega_i`:

.. math::

   \frac{d}{dt} \int_{\Omega_i} \mathbf{U} \, dV
   + \oint_{\partial \Omega_i} \mathbf{F}(\mathbf{U}) \cdot \hat{n} \, dS = 0

This yields an evolution equation for the cell-averaged state vector
:math:`\bar{\mathbf{U}}_i`:

.. math::

   \frac{d \bar{\mathbf{U}}_i}{dt}
   = -\frac{1}{|\Omega_i|} \sum_{\text{faces}} \mathbf{F}(U_L, U_R) \cdot \hat{n} \, A_\text{face}

At each cell interface, a left state :math:`U_L` and a right state :math:`U_R` are reconstructed
from the neighboring cell averages. A Riemann solver then computes the numerical flux
:math:`\mathbf{F}(U_L, U_R)` across the interface. Choosing reconstruction method and Riemann
solver independently allows trading accuracy for robustness.


Reconstruction methods
-----------------------

The reconstruction step recovers high-order approximations of the solution inside each cell from
the cell averages, without introducing spurious oscillations near discontinuities.

**Constant** (1st order)
    Piecewise-constant reconstruction: each cell holds a single value with no variation across it.
    Maximally diffusive but unconditionally robust. Useful for debugging or as a fallback for
    highly irregular data.

**Linear + MinMod** (2nd order)
    A linear slope is fitted inside each cell and then limited by the MinMod limiter, which clips
    the slope to zero whenever neighboring slopes disagree in sign. This prevents spurious
    oscillations near shocks and contact discontinuities while recovering second-order accuracy in
    smooth regions.

**WENO3** (3rd order)
    Weighted Essentially Non-Oscillatory reconstruction of order 3. Nonlinear weights combine
    several candidate stencils so that the scheme automatically degrades to first order near shocks
    (where the smoother stencils receive near-zero weight) while maintaining third-order accuracy
    elsewhere.

**WENOZ** (improved WENO weights)
    An improved weight formula that corrects an accuracy deficiency of WENO3 at smooth extrema.
    WENOZ is less dissipative than WENO3 and generally more accurate on smooth problems, at
    negligible additional cost.

.. tip::

   WENOZ is the recommended default reconstruction for most production runs. It provides
   a good balance between accuracy, robustness near shocks, and computational cost.

**MP5** (5th order monotonicity-preserving)
    A fifth-order scheme with explicit monotonicity constraints that suppress oscillations without
    relying solely on nonlinear weights. MP5 delivers the highest accuracy on smooth solutions and
    is the best choice when the flow is well-resolved and nearly free of discontinuities.


Riemann solvers
----------------

Given the left and right reconstructed states at a cell interface, the Riemann solver returns the
numerical flux. PHARE provides two solvers.

**Rusanov** (local Lax-Friedrichs)
    Uses only the maximum signal speed :math:`S_\max` at the interface:

    .. math::

       \mathbf{F}_\mathrm{Rus} = \frac{1}{2}\left[
         \mathbf{F}(U_L) + \mathbf{F}(U_R) - S_\max (U_R - U_L)
       \right]

    Extremely robust and guaranteed to produce a valid state. The downside is significant
    numerical diffusion. A good starting point when setting up a new problem.

**HLLD**
    Resolves all seven MHD wave families: two fast magnetosonic, two slow magnetosonic, two
    Alfvén, and one entropy wave. By introducing four intermediate states the HLLD solver captures
    contact and rotational discontinuities with far less numerical diffusion than Rusanov.

.. tip::

   HLLD is the recommended Riemann solver for production runs. Use Rusanov only when the problem
   is highly under-resolved or when debugging reconstruction issues.


Time stepping schemes
----------------------

PHARE offers explicit time integrators of increasing order. Higher-order schemes require more
function evaluations per step but allow a larger effective time step without sacrificing stability
or accuracy.

**Euler** (1st order)
    A single forward evaluation per step. Provided for debugging only; the large truncation error
    and tight stability constraint make it unsuitable for production.

**TVDRK2** (2nd order TVD Runge-Kutta)
    A two-stage method that is both second-order accurate and total-variation-diminishing, meaning
    it does not amplify oscillations already present in the solution:

    .. math::

       \mathbf{U}^{(1)} &= \mathbf{U}^n + \Delta t\, L(\mathbf{U}^n) \\
       \mathbf{U}^{n+1} &= \frac{1}{2}\mathbf{U}^n + \frac{1}{2}\left[
         \mathbf{U}^{(1)} + \Delta t\, L(\mathbf{U}^{(1)}) \right]

    where :math:`L` denotes the spatial operator (flux divergence). A solid choice when second-order
    spatial reconstruction (MinMod) is used.

**SSPRK4_5** (4th order, 5 stages)
    A strong stability preserving Runge-Kutta method with four effective orders of accuracy and
    five stages. The additional stages increase the SSP Courant number relative to standard
    four-stage methods, allowing a larger time step for the same stability guarantee.

.. tip::

   SSPRK4_5 is the recommended time integrator when WENOZ or MP5 reconstruction is used. Pairing
   a high-order spatial scheme with a low-order time integrator wastes accuracy and is rarely
   cost-effective.


Dissipative terms
------------------

As in the hybrid PIC formalism (see :doc:`hybridpic`), purely ideal evolution causes current sheets
to sharpen to the grid scale. PHARE supports the same two dissipative mechanisms for the MHD
induction equation.

**Resistivity** :math:`\eta` adds Ohmic diffusion of magnetic flux:

.. math::

   \frac{\partial \mathbf{B}}{\partial t}
   = \nabla \times (\mathbf{v} \times \mathbf{B}) + \eta \nabla^2 \mathbf{B}

**Hyper-resistivity** :math:`\nu` provides scale-selective dissipation that is dominant only at
the smallest resolved scales, leaving large-scale structures intact:

.. math::

   \frac{\partial \mathbf{B}}{\partial t}
   = \nabla \times (\mathbf{v} \times \mathbf{B}) - \nu \nabla^4 \mathbf{B}

These parameters are set in the simulation input (see :doc:`../usage/simulation_inputs`).


Hall MHD
---------

When the Hall term is enabled (``hall=True``), the induction equation is extended to account for
the decoupling of electron and ion flows at sub-ion scales:

.. math::

   \frac{\partial \mathbf{B}}{\partial t}
   = \nabla \times \left( \mathbf{v} \times \mathbf{B}
     - \frac{\mathbf{J} \times \mathbf{B}}{ne} \right)

The additional :math:`\mathbf{J} \times \mathbf{B}` term — the Hall term — introduces whistler
wave dispersion: high-frequency electromagnetic waves whose phase speed increases with wavenumber.
This makes the Hall MHD equations dispersive and tightens the CFL constraint compared to ideal MHD.
Hall MHD is therefore most useful when the domain of interest extends to scales approaching the
ion inertial length :math:`d_i = c/\omega_{pi}`, but a fully kinetic treatment (Hybrid PIC) is
not yet necessary or affordable.


Multi-model coupling
---------------------

PHARE's AMR framework allows different physics solvers to run on different levels of the mesh
hierarchy simultaneously. The parameter ``max_mhd_level`` sets the highest AMR level on which
the MHD solver is used. Levels above that threshold are handled by the Hybrid PIC solver.

This enables a class of simulations in which:

- **Coarse levels** carry the large-scale MHD context (global magnetosphere, solar wind
  conditions), at low computational cost.
- **Fine levels** resolve kinetic physics in regions of interest — for example, a magnetic
  reconnection site or a collisionless shock — with full particle dynamics.

The two solvers exchange boundary conditions across the MHD–PIC interface automatically through
PHARE's messenger infrastructure (see :doc:`amr`).

.. tip::

   Multi-model coupling is particularly powerful for global magnetospheric simulations: the
   MHD solver handles the magnetosheath and tail lobes while PIC resolves the diffusion region of
   a reconnection event embedded within the global structure.
