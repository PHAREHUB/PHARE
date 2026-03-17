
Tutorial: 2D Orszag-Tang Vortex
================================

The Orszag-Tang vortex is a canonical 2D MHD benchmark. Smooth sinusoidal
initial conditions evolve into a network of interacting shocks, current sheets,
and vortices. Because the exact nonlinear solution is not known analytically,
codes are compared against each other and against high-resolution reference runs.
It is the go-to problem for checking that a multidimensional MHD solver handles
shock interactions, current-sheet steepening, and the divergence-free constraint
on **B** correctly.

.. contents:: Contents
   :local:
   :depth: 2


The physics
-----------

The problem lives on a doubly periodic square domain of side :math:`L`. The
initial velocity and magnetic fields are both purely sinusoidal:

.. math::

   v_x = -\sin\!\left(\frac{2\pi y}{L}\right), \qquad
   v_y =  \sin\!\left(\frac{2\pi x}{L}\right)

.. math::

   B_x = -B_0 \sin\!\left(\frac{2\pi y}{L}\right), \qquad
   B_y =  B_0 \sin\!\left(\frac{4\pi x}{L}\right)

with :math:`B_0 = 1/\!\sqrt{4\pi}`. The velocity field is an incompressible
vortex; the magnetic field has a different spatial frequency in *x* and *y*,
so the two patterns are incommensurate and immediately drive nonlinear wave
interactions.

Density and pressure are uniform at :math:`t = 0`:

.. math::

   \rho_0 = \frac{25}{36\pi}, \qquad p_0 = \frac{5}{12\pi}

These values set the plasma beta and Alfven Mach number to order unity, so
neither the kinetic nor magnetic pressure dominates.

The nonlinear evolution follows a well-known sequence. By :math:`t \approx 0.5`
(in the natural units of the problem) the smooth initial fields have steepened
into a pair of interacting fast shocks. By :math:`t \approx 1` a coherent
current sheet forms near the center of the domain and current-sheet
reconnection begins. At late times the solution develops a turbulent-like
cascade of structures at smaller and smaller scales, making the problem a
sensitive test of numerical dissipation and stability.

This tutorial runs PHARE's MHD solver with Linear reconstruction, MinMod
limiting, and a Rusanov Riemann solver on a 256 x 256 grid with Hall MHD and
hyper-resistivity for numerical stability.

For the theoretical background on the MHD model implemented in PHARE, see
:doc:`../theory/mhd`.


Setting up the simulation
--------------------------

Create a file ``orszag_tang_2d.py`` with the content below. Each parameter
choice is explained in the inline comments.

.. code-block:: python

    #!/usr/bin/env python3

    import numpy as np
    import pyphare.pharein as ph
    from pyphare.simulator.simulator import Simulator, startMPI

    ph.NO_GUI()

    final_time    = 25.6
    time_step     = 0.0005
    diag_dir      = "phare_outputs/orszag_tang"

    time_step_nbr = int(final_time / time_step)

    # Write diagnostics every 1000 steps (every 0.5 time units).
    dumpfrequency = 1000
    dt_dump       = dumpfrequency * time_step
    timestamps    = (
        dt_dump * np.arange(int(final_time / dt_dump) + 1)
    )


    def config():
        # ------------------------------------------------------------------
        # Domain: 256 x 256 cells of size 0.2 => physical side L = 51.2.
        # The doubly periodic domain requires cells * dl to be the same in
        # both directions for the standard Orszag-Tang initial conditions.
        # ------------------------------------------------------------------
        cells = (256, 256)
        dl    = (0.2, 0.2)

        sim = ph.Simulation(
            smallest_patch_size=15,
            time_step_nbr=time_step_nbr,
            time_step=time_step,
            cells=cells,
            dl=dl,
            interp_order=2,
            # ------------------------------------------------------------------
            # AMR: single MHD level for this tutorial.  Tagging-based refinement
            # is still registered so the code path is exercised, but with
            # max_nbr_levels=1 no additional levels are created.
            # ------------------------------------------------------------------
            refinement="tagging",
            max_mhd_level=1,
            max_nbr_levels=1,
            nesting_buffer=1,
            # ------------------------------------------------------------------
            # Fluid MHD model.
            # ------------------------------------------------------------------
            model_options=["MHDModel"],
            gamma=5.0 / 3.0,
            # ------------------------------------------------------------------
            # Reconstruction, limiter, and Riemann solver.
            # Linear reconstruction with MinMod limiting gives second-order
            # spatial accuracy while remaining TVD (no new extrema are created).
            # Rusanov is the simplest approximate Riemann solver; it is robust
            # across all wave families at the cost of moderate numerical
            # diffusion.
            # ------------------------------------------------------------------
            reconstruction="Linear",
            limiter="MinMod",
            riemann="Rusanov",
            # ------------------------------------------------------------------
            # Time integration: TVDRK2 is a second-order strong-stability-
            # preserving Runge-Kutta scheme that keeps the TVD property of the
            # reconstruction step.
            # ------------------------------------------------------------------
            mhd_timestepper="TVDRK2",
            # ------------------------------------------------------------------
            # Dissipation.
            # Resistivity (eta, resistivity) is zero for this ideal run.
            # A small kinematic hyper-resistivity (nu=0.02) damps grid-scale
            # noise introduced by the shock interactions without noticeable
            # impact on the macroscopic flow structure.
            # ------------------------------------------------------------------
            eta=0.0,
            resistivity=0.0,
            hyper_resistivity=0.0,
            nu=0.02,
            hyper_mode="spatial",
            # Hall MHD term: couples ion and electron momentum at sub-ion scales.
            hall=True,
            res=False,
            hyper_res=True,
            diag_options={
                "format": "phareh5",
                "options": {"dir": diag_dir, "mode": "overwrite"},
            },
            strict=True,
        )

        # ------------------------------------------------------------------
        # Amplitude of the magnetic field.
        # B0 = 1/sqrt(4 pi) sets the Alfven speed to approximately 0.28 in
        # simulation units and keeps the problem in the standard Orszag-Tang
        # parameter regime.
        # ------------------------------------------------------------------
        B0 = 1.0 / np.sqrt(4.0 * np.pi)

        # Physical domain lengths are available after the Simulation object is
        # built (they equal cells * dl).
        Lx = sim.simulation_domain()[0]   # 51.2
        Ly = sim.simulation_domain()[1]   # 51.2

        # ------------------------------------------------------------------
        # Initial conditions.
        # ------------------------------------------------------------------
        def density(x, y):
            return 25.0 / (36.0 * np.pi)

        def vx(x, y):
            return -np.sin(2.0 * np.pi * y / Ly)

        def vy(x, y):
            return np.sin(2.0 * np.pi * x / Lx)

        def vz(x, y):
            return 0.0

        def bx(x, y):
            return -B0 * np.sin(2.0 * np.pi * y / Ly)

        def by(x, y):
            return B0 * np.sin(4.0 * np.pi * x / Lx)

        def bz(x, y):
            return 0.0

        def p(x, y):
            return 5.0 / (12.0 * np.pi)

        ph.MHDModel(density=density, vx=vx, vy=vy, vz=vz,
                    bx=bx, by=by, bz=bz, p=p)

        # ------------------------------------------------------------------
        # Diagnostics: electromagnetic field B and MHD fluid quantities.
        # All quantities are written at every dump timestamp.
        # ------------------------------------------------------------------
        ph.ElectromagDiagnostics(quantity="B", write_timestamps=timestamps)

        for quantity in ["rho", "V", "P"]:
            ph.MHDDiagnostics(quantity=quantity, write_timestamps=timestamps)

        return sim


    if __name__ == "__main__":
        startMPI()
        Simulator(config()).run()

For a complete reference of all accepted parameters, see
:doc:`../usage/simulation_inputs`.


Key parameter choices explained
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Parameter
     - Value
     - Rationale
   * - ``cells``
     - ``(256, 256)``
     - Standard resolution for the Orszag-Tang benchmark; resolves the main
       shock structures while remaining affordable on a workstation.
   * - ``dl``
     - ``(0.2, 0.2)``
     - Cell size in simulation units; gives a physical domain of 51.2 x 51.2.
   * - ``gamma``
     - ``5/3``
     - Monatomic ideal gas; required by the reference solution.
   * - ``reconstruction``
     - ``"Linear"``
     - Second-order piecewise-linear reconstruction of interface states.
   * - ``limiter``
     - ``"MinMod"``
     - Slope limiter that enforces TVD: prevents new extrema from forming
       and keeps the scheme stable near shocks.
   * - ``riemann``
     - ``"Rusanov"``
     - Local Lax-Friedrichs Riemann solver. Robust for all MHD wave families;
       introduces moderate numerical diffusion.
   * - ``mhd_timestepper``
     - ``"TVDRK2"``
     - Second-order strong-stability-preserving Runge-Kutta. Preserves the TVD
       property of the spatial reconstruction step.
   * - ``nu``
     - ``0.02``
     - Hyper-resistivity coefficient. Damps grid-scale oscillations generated
       by shock interactions without affecting large-scale dynamics.
   * - ``hall``
     - ``True``
     - Activates the Hall term in the induction equation. Relevant at the
       sub-ion scales that arise in the late-time current sheets.


Running the simulation
-----------------------

From the directory containing ``orszag_tang_2d.py``, run:

.. code-block:: bash

    mpirun -np 4 python3 -Ou orszag_tang_2d.py

The ``-O`` flag disables Python assertion overhead and ``-u`` forces unbuffered
output so progress messages appear in real time.

The 256 x 256 grid runs 51 200 time steps in total. On a four-core workstation
expect a wall-clock time of roughly 30-60 minutes depending on hardware.
Diagnostic files accumulate in ``phare_outputs/orszag_tang/`` as HDF5 archives
(``EM_B.h5``, ``mhd_rho.h5``, ``mhd_V.h5``, ``mhd_P.h5``).

For a quick sanity check, reduce ``cells`` to ``(64, 64)`` and increase
``time_step`` to ``0.002``; the run completes in a few minutes and the main
shock structures are still visible.


Analyzing results
-----------------

Use ``pharesee`` to load the HDF5 output and plot 2D colormaps of the density
at several timestamps. Density is one of the clearest diagnostics because it
highlights both the shock fronts (density jumps) and the low-density regions
carved out by the vortex.

.. code-block:: python

    import numpy as np
    import matplotlib.pyplot as plt
    from pyphare.pharesee.run import Run

    diag_dir = "phare_outputs/orszag_tang"
    run = Run(diag_dir)

    # Pick three snapshots: early (smooth), intermediate (shock formation),
    # and late (complex turbulent-like state).
    selected_times = [0.5, 5.0, 25.0]

    fig, axes = plt.subplots(1, 3, figsize=(15, 5), constrained_layout=True)

    for ax, t in zip(axes, selected_times):
        rho = run.GetMHDrho(t, merged=True)
        xi, yi = rho["rho"][1]
        rho_vals = rho["rho"][0](xi, yi)

        pcm = ax.pcolormesh(xi, yi, rho_vals.T, cmap="inferno", shading="auto")
        ax.set_title(f"$\\rho$,  t = {t:.1f}")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_aspect("equal")
        fig.colorbar(pcm, ax=ax, label=r"$\rho$")

    fig.suptitle("Orszag-Tang vortex — density evolution")
    plt.show()

``run.GetMHDrho(t, merged=True)`` returns a dictionary with a single key
``"rho"``. The value is a tuple ``(interpolant, [x_coords, y_coords])``.
Calling ``interpolant(xi, yi)`` evaluates the density on the coordinate
arrays ``xi`` and ``yi``.

Magnetic field components are accessed in the same way through
``run.GetB(t, merged=True)``, which returns keys ``"Bx"``, ``"By"``, and
``"Bz"``.

For the full pharesee data API, see :doc:`../pharesee/get_data`.


Plotting divergence of B
~~~~~~~~~~~~~~~~~~~~~~~~

The :math:`\nabla \cdot \mathbf{B} = 0` constraint is maintained numerically
by the constrained transport scheme. Monitoring its violation is a standard
quality check for MHD simulations:

.. code-block:: python

    import matplotlib.pyplot as plt
    from pyphare.pharesee.run import Run

    run = Run("phare_outputs/orszag_tang")
    t   = 25.0

    divb = run.GetDivB(t)
    divb.plot(
        filename="divb_t25.png",
        plot_patches=True,
        vmin=-1e-11,
        vmax=1e-11,
    )

``GetDivB`` computes the cell-centered divergence from the face-centered field
values stored by the constrained-transport integrator. Values should remain
at machine-epsilon level throughout the run.


Expected output
---------------

The density field passes through three qualitatively distinct phases:

1. **Early time** (:math:`t \lesssim 1`): the sinusoidal density profile
   is still recognizable but the regions of strong velocity shear steepen
   and the first fast shocks appear as sharp density ridges.

2. **Intermediate time** (:math:`t \sim 5{-}15`): interacting shocks
   from opposite sides of the domain collide near the domain center,
   producing a concentrated current sheet and a pronounced density enhancement.
   Circular rarefaction regions form around the original vortex centers.

3. **Late time** (:math:`t \sim 25`): the current sheet has fragmented
   into multiple secondary structures. The density map shows a complex,
   pseudo-turbulent pattern with structures at many length scales.

The magnetic pressure :math:`B^2/2` traces similar features offset by
the Alfven travel time. The pressure :math:`P` shows the strongest signatures
at shock intersections.

.. note::

   This run uses a single AMR level (``max_nbr_levels=1``), so no adaptive
   refinement is applied. The late-time current sheets are therefore resolved
   only at the coarse 256 x 256 resolution. Enabling a second AMR level
   (``max_nbr_levels=2``) allows the code to track the thinnest current
   structures with higher resolution cells while keeping the large-scale grid
   affordable.


Next steps
----------

- Increase ``max_nbr_levels`` to 2 to let AMR follow the thinnest current
  sheets automatically. Compare the refined and unrefined density maps.
- Switch ``riemann`` from ``"Rusanov"`` to ``"hlld"`` for a less diffusive
  flux computation and sharper shock transitions.
- Set ``hall=False`` and compare the late-time current-sheet morphology to
  assess the role of the Hall term at the sub-ion scales.
- Add a second output quantity: velocity magnitude
  :math:`|\mathbf{v}| = \sqrt{v_x^2 + v_y^2}` to visualize the vortex
  structure alongside the density shocks.

See also:

- :doc:`../usage/simulation_inputs` — full reference for all configuration blocks
- :doc:`../pharesee/get_data` — reading and manipulating simulation output
- :doc:`../theory/mhd` — theoretical background of the MHD model
