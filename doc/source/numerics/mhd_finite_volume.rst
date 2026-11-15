=========================
MHD finite-volume scheme
=========================

**Goal:** understand how PHARE's MHD solver computes a step, and choose its numerical
options: reconstruction, limiter, Riemann solver, Hall term and dissipation.

The equations solved are on :doc:`/theory/mhd`. How to write an MHD input script is on
:doc:`/user_guide/mhd_inputs`.

.. contents:: On this page
   :local:
   :depth: 1

The details of each part have their own page:
:doc:`mhd_reconstruction` (values at the faces, limiters),
:doc:`mhd_riemann_solvers` (the face flux and the electric field on the edges),
:doc:`mhd_time_integration` (Runge–Kutta integrators and time step) and
:doc:`/validation/mhd` (which tests exist and what was measured).

Minimal example
---------------

The scheme is chosen in the ``Simulation`` block. This is a fragment; a complete,
runnable MHD script is on :doc:`/user_guide/mhd_inputs` (the 1D shock-tube test uses
exactly these choices):

.. code-block:: python

   sim = ph.Simulation(
       model_options=["MHDModel"],
       reconstruction="Linear",   # how cell averages are extended to the faces
       limiter="VanLeer",         # only used by "Linear"
       riemann="Rusanov",         # flux at each face
       mhd_timestepper="TVDRK2",  # time integrator, see numerics/time_integration
       hall=False,
       gamma=5.0 / 3.0,
       eta=0.0,                   # resistivity
       nu=0.0,                    # hyper-resistivity (needs hall=True)
       # ... grid, time and diagnostics options
   )

What it does
------------

PHARE's MHD solver is a :term:`Godunov method` (a :term:`finite volume` scheme). The fluid variables (density :math:`\rho`, momentum :math:`\rho\mathbf{v}`, total energy
:math:`E_{tot}`) are *cell averages*, stored at cell centres. The magnetic field is stored
on the cell faces, one component per face, and the electric field on the cell edges
(:doc:`grid`). One stage of the scheme does:

#. **Reconstruction**: from the cell averages, two values at each face, one from each side.
#. :term:`Riemann solver`: from that jump, the flux of mass, momentum, energy and B through
   the face.
#. **Finite-volume update**: :math:`U_i \mathrel{-}= \frac{\Delta t}{\Delta x}(F_{i+1/2} - F_{i-1/2})`
   in each direction.
#. :term:`Constrained transport <constrained transport>`: E on the cell edges from the
   Riemann solvers of the faces meeting there (an *upwind* average), then
   :math:`\partial_t\mathbf{B} = -\nabla\times\mathbf{E}` on the faces.

.. plot::
   :include-source: False
   :caption: Why constrained transport keeps ∇·B (more on :doc:`div_b`). In 2D, Ez sits on the four corners of a
      cell; each corner value changes the two faces that meet there, with opposite signs in
      the cell's balance, so the net flux of B out of the cell does not change.

   from phare_draw import *

   fig, ax = strip(1, phone=True)
   a = ax[0]
   grid(a, (0, 0, 1, 1), color=C["grey"], lw=1.2)
   for (x, y) in [(0, 0), (1, 0), (0, 1), (1, 1)]:
       mark(a, "primal", x, y, color=C["E"])
   text(a, -0.12, 1.1, "Ez", color=C["E"], bold=True)
   for x in (0, 1):
       arrow(a, (x - 0.15, 0.5), (x + 0.15, 0.5), color=C["B"], lw=2.2, ms=12)
       arrow(a, (x, 0.08), (x, 0.35), color=C["grey"], lw=1.2, ms=9)
       arrow(a, (x, 0.92), (x, 0.65), color=C["grey"], lw=1.2, ms=9)
   for y in (0, 1):
       arrow(a, (0.5, y - 0.15), (0.5, y + 0.15), color=C["B"], lw=2.2, ms=12)
       arrow(a, (0.08, y), (0.35, y), color=C["grey"], lw=1.2, ms=9)
       arrow(a, (0.92, y), (0.65, y), color=C["grey"], lw=1.2, ms=9)
   text(a, 1.2, 0.5, "Bx", color=C["B"]); text(a, 0.5, 1.25, "By", color=C["B"])
   text(a, 0.5, 0.5, "∇·B\nunchanged", color=C["grey"], fs=11)
   a.update_datalim([(-0.4, -0.2), (1.45, 1.4)])
   finish(fig, ax)

.. plot::
   :include-source: False
   :caption: Reconstruction in 1D. Top: constant (first order), each cell is flat and the
      jump at every face is the full difference between neighbours. Bottom: linear with a
      slope limiter (second order), each cell gets a slope, so the left and right values
      at a face are much closer on smooth parts; the limiter flattens the slope next to
      the discontinuity to avoid creating new extrema.

   import numpy as np
   import matplotlib.pyplot as plt

   u = np.array([1.0, 1.0, 0.95, 0.8, 0.55, 0.3, 0.15, 0.1, 0.1, 0.1])
   u[6:] = [0.9, 0.9, 0.9, 0.9]  # a jump
   n = len(u)
   fig, axes = plt.subplots(2, 1, figsize=(3.4, 4.2), sharex=True)
   for ax, title in zip(axes, ["Constant", "Linear + van Leer limiter"]):
       for i in range(n):
           if title == "Constant":
               s = 0.0
           else:
               dl = u[i] - u[i - 1] if i > 0 else 0.0
               dr = u[i + 1] - u[i] if i < n - 1 else 0.0
               s = 2 * dl * dr / (dl + dr) if dl * dr > 0 else 0.0
           ax.plot([i, i + 1], [u[i] - s / 2, u[i] + s / 2], color="#2a7ab9", lw=2.5)
           ax.plot(i + 0.5, u[i], "o", color="0.3", ms=5)
       for i in range(n + 1):
           ax.axvline(i, color="0.85", lw=0.8, zorder=0)
       ax.set_title(title, fontsize=12)
       ax.set_yticks([])
       ax.tick_params(labelsize=11)
   axes[1].set_xlabel("cell faces", fontsize=12)
   fig.tight_layout()

One stage on a row of cells
---------------------------

.. plot::
   :include-source: False
   :alt: One stage of the scheme in 1D: cell values, primitive values, left and right value at a face, flux, update of a cell.
   :caption: One stage in 1D. 1: the stored cell values (density, momentum, total energy; B on the faces). 2: density, velocity and pressure, cell by cell. 3: the reconstruction gives a left and a right value at each face. 4: the Riemann solver gives the flux through the face. 5: each cell is updated from the two fluxes around it.

   from phare_draw import *

   fig, ax = strip(5, ["Cell values", "Primitive values", "Left and right at a face", "Flux through a face", "Update of a cell"], phone=True, height=1.5)
   n = 5
   for k, a in enumerate(ax):
       grid(a, (0, 0, n, 1), color=C["grey"], lw=1.2)
       for i in range(n):
           mark(a, "primal", i + 0.5, 0.5, color=C["ion"] if k < 2 else C["light"], size=7)
   x = 3.0
   mark(ax[2], "primal", x - 0.25, 0.5, color=C["L1"], size=8); text(ax[2], x - 0.55, 1.2, "left", color=C["L1"], fs=9)
   mark(ax[2], "primal", x + 0.25, 0.5, color=C["L2"], size=8); text(ax[2], x + 0.55, 1.2, "right", color=C["L2"], fs=9)
   arrow(ax[3], (x - 0.6, 0.5), (x + 0.6, 0.5), color=C["B"], lw=2.4, ms=12); text(ax[3], x, 1.25, "flux F", color=C["B"], fs=9)
   arrow(ax[4], (1.9, 0.5), (2.4, 0.5), color=C["B"], lw=2.2, ms=10); arrow(ax[4], (3.4, 0.5), (3.9, 0.5), color=C["B"], lw=2.2, ms=10)
   mark(ax[4], "primal", 3.0 + 0.0, 0.5, color=C["ion"], size=9)
   text(ax[4], 2.9, 1.3, "new = old − Δt/Δx (F_right − F_left)", fs=8)
   finish(fig, ax)

Options
-------

* ``reconstruction``: ``"Constant"``, ``"Linear"``, ``"WENO3"``, ``"WENOZ"``, ``"MP5"``.
  WENO: weighted essentially non-oscillatory; MP5: fifth-order monotonicity-preserving. The
  scheme is second order at best whatever the reconstruction (see `Order of the scheme`_): a
  wider reconstruction lowers the numerical diffusion, not the order. Ghost layers: 2 for
  ``Constant``, 4 for the others. No default: must be given. Details:
  :doc:`mhd_reconstruction`.
* ``limiter``: ``"None"``, ``"VanLeer"``, ``"MinMod"``. Only ``Linear`` uses a limiter, and it
  needs one: ``VanLeer`` or ``MinMod``. All other reconstructions need ``"None"``. No default:
  must be given. ``MinMod`` is more diffusive than ``VanLeer``.
* ``riemann``: ``"Rusanov"``, ``"HLL"``, ``"HLLD"``. ``Rusanov`` uses one wave speed (the
  fastest) and is the most diffusive. ``HLL`` (a two-wave approximate solver) uses the
  fastest left- and right-going speeds. ``HLLD`` (HLL with intermediate discontinuities) also
  resolves the Alfvén and contact waves, so it is the least diffusive. The wave speeds are
  fast-magnetosonic only, without the whistler speed. No default: must be given. Details:
  :doc:`mhd_riemann_solvers`.
* ``hall``: ``True`` or ``False``. Adds the Hall term :math:`\mathbf{J}\times\mathbf{B}/\rho` to
  the electric field (Hall MHD; :math:`\rho = n` in PHARE units). Its whistler waves make the
  time step scale like :math:`\Delta x^2`.
* ``gamma``: number > 1. Adiabatic index. Default 5/3.
* ``eta``: number ≥ 0. Resistivity :math:`\eta`. Default 0.
* ``nu``: number ≥ 0. Hyper-resistivity :math:`\nu`. Default 0. Requires ``hall=True``.
* ``hyper_mode``: ``"constant"`` or ``"spatial"``. ``"spatial"`` uses
  :math:`\nu\,\Delta x^2\,(|\mathbf{B}|/\rho + 1)` with :math:`\Delta x` the smallest cell
  size, so the dissipation follows the grid. Default ``"constant"``.

``reconstruction``, ``limiter`` and ``riemann`` are written exactly as above
(case-sensitive). ``mhd_timestepper`` is described on :doc:`mhd_time_integration`.

Order of the scheme
-------------------

**PHARE's MHD scheme is second order at best, in 1D, 2D and 3D.** The reason, in one line:
the stored cell values are used as if they were point values, and each face gets one flux
value, so the scheme cannot be better than second order whatever the reconstruction
(:src:`src/amr/solvers/time_integrator/compute_fluxes.hpp#L66`,
:src:`src/core/numerics/godunov_fluxes/godunov_fluxes.hpp#L106-L146`). This is an argument from the
code. The only measurement above one dimension is the 2D AMR convergence tests, which accept
orders between 1.70 and 2.25 (:doc:`/validation/mhd`).

* The reconstructions are chosen for **less numerical diffusion and sharper features** on a
  given grid. Their own orders (1, 2, 3, 5, 5) are properties of the reconstructions,
  visible on a small linear wave (:doc:`mhd_reconstruction`); they are not the order of the
  scheme.
* The integrators have orders 1 to 4 in time (:doc:`mhd_time_integration`). The result is no
  better than the lowest order of its parts.

.. admonition:: In preparation
   :class: note

   A fourth-order version of the MHD scheme is in development.

How to choose
-------------

Only what the tests and measurements of this documentation support:

* **Shocks and discontinuities, as in the tests of the repository.** The shock-tube and
  Orszag–Tang tests use ``Linear`` + ``VanLeer`` + ``Rusanov`` + ``TVDRK2``
  (:src:`tests/functional/mhd_shock/common.py#L77-L80`,
  :src:`tests/functional/mhd_orszagtang/orszag_tang.py#L59-L62`). On the shock tube,
  ``Linear`` + ``VanLeer`` smears the contact over 7 cells, ``Constant`` + ``Rusanov`` over
  45 (:doc:`mhd_riemann_solvers`).
* **Less diffusion for the same reconstruction.** ``HLLD`` smears the contact about half as
  much as ``Rusanov`` on that test; it is built in 1D with ``Constant`` only in the default
  list (:ref:`below <which-combinations>`).
* **Smooth flows.** ``WENOZ`` or ``MP5`` give the smallest left–right difference on a
  smooth profile (:doc:`mhd_reconstruction`); pair them with ``TVDRK3`` or ``SSPRK4_5``
  (a recommendation of the literature for WENO-Z; PHARE does not enforce it; the tests that use ``WENOZ``
  do so).
* **Cost.** Ghost layers: 2 for ``Constant``, 4 for the others. Flux computations per
  step: 1, 2, 3, 5 for ``Euler``, ``TVDRK2``, ``TVDRK3``, ``SSPRK4_5``.
* **Hall.** With ``hall=True`` the step must follow the whistler limit, :math:`\Delta t\propto\Delta x^2`
  (:doc:`mhd_time_integration`).

Extra terms
-----------

* **Hall** (``hall=True``): the field :math:`\mathbf J\times\mathbf B/\rho` is added to the
  flux of B and of energy, and to the electric field on the edges
  (:src:`src/core/numerics/MHD_equations/MHD_equations.hpp#L113-L144`).
* **Resistivity** ``eta``: :math:`\mathbf E_\mathrm{diss}=\eta\mathbf J`. **Hyper-resistivity** ``nu``:
  :math:`-\nu\nabla^2\mathbf J`, or with ``hyper_mode="spatial"`` the coefficient
  :math:`\nu\,\Delta x^2(|\mathbf B|/\rho+1)`
  (:src:`src/core/numerics/constrained_transport/dissipative_electric_field.hpp#L63-L113`).
  They add to the electric field and to the fluxes of B and energy.
* There is no density or pressure floor: the line that would apply one is commented out
  (:src:`src/core/numerics/primite_conservative_converter/to_primitive_converter.hpp#L26-L28`).

.. _which-combinations:

Which combinations exist
------------------------

``reconstruction``, ``limiter``, ``riemann``, ``hall`` and the dimension are fixed when
PHARE is compiled: each combination is a separate compiled module. The list built by
default is in :src:`res/sim/all.txt`; among others it contains, in 1D, 2D and 3D,
``Linear`` + ``VanLeer`` + ``Rusanov`` and ``WENOZ`` + ``None`` + ``Rusanov``. ``HLL`` and
``HLLD`` are built in 1D only, with ``Constant``.

If you ask for a combination that was not built, the run stops at start-up with a
message like::

   No module named 'pybindlibs.cpp_2_MP5_None_HLLD_false'.
   This run needs the build permutation
       2,MP5,None,HLLD,false
   Add it to a res/sim permutation file (res/sim/all.txt) and rebuild.

followed by the list of combinations that are built. Add the printed line to the file
(or to your own file, passed to CMake with ``-DPHARE_PERMUTATIONS=<file>``) and rebuild.

Common mistakes
---------------

* **A limiter with a non-linear reconstruction, or no limiter with** ``Linear`` — the
  combination cannot be compiled — use ``"None"`` except with ``Linear``.
* **Wrong case, for example** ``"weno3"`` — "No module named ..." at start-up — use the
  exact names above.
* **Setting** ``hyper_resistivity`` **in an MHD run** — no effect on MHD — that option is the
  hybrid one; MHD uses ``nu``.
* ``nu > 0`` **without** ``hall=True`` — the run stops with
  ``hyper-resistivity (nu > 0) requires the Hall term to be enabled``.
* **Expecting a density or pressure floor** — there is none: a negative pressure is not
  corrected, and in a release build the run may continue with NaN values, so check your
  diagnostics. Reduce the time step, or use a more diffusive reconstruction, limiter or
  Riemann solver near strong shocks.

.. seealso:: How this works inside: :doc:`/dev/mhd_solver`.
