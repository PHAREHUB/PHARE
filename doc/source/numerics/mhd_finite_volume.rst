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

.. mermaid::

   flowchart TD
     A["cell averages ρ, ρv, Etot<br/>(+ B on faces)"] --> R["reconstruct left/right values<br/>at every face"]
     R --> F["Riemann solver<br/>flux through every face"]
     F --> U["update the cell averages<br/>with the flux differences"]
     F --> CT["constrained transport<br/>E on every edge"]
     CT --> B["Faraday: update B on faces<br/>with the curl of E"]

1. **Reconstruction.** Inside each cell, the solution is rebuilt from the averages of the
   cell and its neighbours, which gives two values at each face: one from the left cell,
   one from the right cell.
2. :term:`Riemann solver`. At each face, the jump between the two values is treated as a small
   Riemann problem; its approximate solution gives the flux of mass, momentum, energy and
   magnetic field through the face.
3. **Finite-volume update.** Each cell average changes by what flows in minus what flows
   out: :math:`U_i^{new} = U_i - \frac{\Delta t}{\Delta x}(F_{i+1/2} - F_{i-1/2})`, in each
   direction.
4. :term:`Constrained transport <constrained transport>`. B is not updated with its own face fluxes. Instead, the
   electric field is built on the cell edges from the information of the Riemann
   solvers of the faces that meet at that edge (an *upwind* average), and B is advanced
   with Faraday's law, :math:`\partial_t\mathbf{B} = -\nabla\times\mathbf{E}`. Since the
   discrete divergence of a discrete curl is exactly zero on this grid (to round-off), the update never
   changes :math:`\nabla\cdot\mathbf{B}`: no divergence cleaning is needed.

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
   fig, axes = plt.subplots(2, 1, figsize=(5, 5), sharex=True)
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
       ax.set_title(title, fontsize=15)
       ax.set_yticks([])
       ax.tick_params(labelsize=14)
   axes[1].set_xlabel("cell faces", fontsize=15)
   fig.tight_layout()

Options
-------

.. list-table::
   :header-rows: 1
   :widths: 18 40 42

   * - Option
     - Values
     - Notes
   * - ``reconstruction``
     - ``"Constant"``, ``"Linear"``, ``"WENO3"``, ``"WENOZ"``, ``"MP5"``
     - WENO: weighted essentially non-oscillatory; MP5: fifth-order
       monotonicity-preserving. In 1D, orders on smooth flows: 1, 2, 3, 5, 5. In 2D and 3D
       the whole scheme is second order (one flux point per face), so higher-order
       reconstructions lower the error, not the order. Wider stencils need more ghosts:
       2 layers for ``Constant``, 4 for the others. No default: must be given.
   * - ``limiter``
     - ``"None"``, ``"VanLeer"``, ``"MinMod"``
     - Only ``Linear`` uses a limiter, and it needs one: ``VanLeer`` or ``MinMod``. All other
       reconstructions need ``"None"``. No default: must be given. ``MinMod`` is more diffusive than ``VanLeer``.
   * - ``riemann``
     - ``"Rusanov"``, ``"HLL"``, ``"HLLD"``
     - ``Rusanov`` uses one wave speed (the fastest) and is the most diffusive. ``HLL``
       (a two-wave approximate solver) uses the fastest left- and right-going speeds.
       ``HLLD`` (HLL with intermediate discontinuities) also resolves the Alfvén and
       contact waves, so it is the least diffusive. The wave speeds are fast-magnetosonic
       only, without the whistler speed. No default: must be given.
   * - ``hall``
     - ``True``, ``False``
     - Adds the Hall term :math:`\mathbf{J}\times\mathbf{B}/\rho` to the electric field
       (Hall MHD; :math:`\rho = n` in PHARE units). Its whistler waves make the time step scale like :math:`\Delta x^2`.
   * - ``gamma``
     - number > 1
     - Adiabatic index. Default 5/3.
   * - ``eta``
     - number ≥ 0
     - Resistivity :math:`\eta`. Default 0.
   * - ``nu``
     - number ≥ 0
     - Hyper-resistivity :math:`\nu`. Default 0. Requires ``hall=True``.
   * - ``hyper_mode``
     - ``"constant"``, ``"spatial"``
     - ``"spatial"`` uses :math:`\nu\,\Delta x^2\,(|\mathbf{B}|/\rho + 1)` with
       :math:`\Delta x` the smallest cell size, so the dissipation follows the grid.
       Default ``"constant"``.

``reconstruction``, ``limiter`` and ``riemann`` are written exactly as above
(case-sensitive). ``mhd_timestepper`` is described on :doc:`time_integration`.

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
