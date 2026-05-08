
Tutorial: 2D Harris Sheet Reconnection
=======================================

This tutorial walks through a 2D magnetic reconnection simulation using the
Harris current sheet equilibrium. It is the standard benchmark for hybrid PIC
codes studying collisionless reconnection and introduces several features absent
from the 1D Alfven wave tutorial: 2D initialization functions, the double
current sheet geometry, hyper-resistivity, and reconnection-specific analysis
tools from ``pharesee``.

.. contents:: Contents
   :local:
   :depth: 2


The physics
-----------

Magnetic reconnection is a fundamental energy-conversion process in space and
astrophysical plasmas. It occurs wherever oppositely directed magnetic field
lines are brought together, as at the dayside magnetopause, in the magnetotail,
and in solar flares. At the reconnection site the field topology changes,
converting magnetic energy into particle kinetic energy through the formation of
fast outflow jets and energetic particle populations.

The Harris current sheet is the textbook equilibrium for studying reconnection.
In its ideal 1D form the magnetic field reverses direction across a thin current
layer of half-width *L*:

.. math::

   B_x(y) = B_0 \tanh\!\left(\frac{y - y_0}{L}\right)

The self-consistent Harris density profile that maintains pressure balance is:

.. math::

   n(y) = \frac{n_0}{\cosh^2\!\!\left(\dfrac{y - y_0}{L}\right)} + n_\mathrm{bg}

where :math:`n_0` is the peak current-sheet density and :math:`n_\mathrm{bg}` is
a uniform background density that keeps particle statistics finite away from the
sheet.

The total pressure :math:`p + B^2/2` must be uniform for the equilibrium to
hold, which ties the local temperature to the local density and field strength:

.. math::

   T(x, y) = \frac{K - B^2(x,y)/2}{n(x,y)}

where *K* is the total pressure constant.

The equilibrium is linearly unstable to the *tearing mode*: a resistive
instability that spontaneously breaks and reconnects field lines, forming a
chain of magnetic islands (plasmoids) at the neutral line :math:`B_x = 0`. In
the nonlinear phase the islands coalesce and reconnection drives strong
out-of-plane current jets (detectable through :math:`J_z`) and characteristic
X-point topologies in the magnetic flux function :math:`A_z`.

In this tutorial a double Harris sheet is used, placing two neutral lines at
:math:`y = 0.3 L_y` and :math:`y = 0.7 L_y`. A small magnetic perturbation
seeds the tearing mode so that reconnection starts within a reasonable
simulation time. For background on the hybrid PIC model see
:doc:`../theory/hybridpic`, and for the AMR strategy see :doc:`../theory/amr`.


Setting up the simulation
--------------------------

Create a file ``harris_2d.py`` with the content below. Each parameter choice
is explained in the inline comments.

.. code-block:: python

    #!/usr/bin/env python3

    import numpy as np
    import pyphare.pharein as ph
    from pyphare.simulator.simulator import Simulator, startMPI

    ph.NO_GUI()


    diag_dir = "phare_outputs/harris_2d"

    # ----------------------------------------------------------------
    # Timestamps at which diagnostics are written.
    # For a production run set these to, e.g.:
    #   np.arange(0, final_time + dt, final_time / 20)
    # ----------------------------------------------------------------
    final_time = 50.0
    time_step  = 0.005
    timestamps = np.arange(0, final_time + time_step, final_time / 20)


    def config():
        # Half-width of each current sheet in ion inertial lengths.
        L = 0.5

        sim = ph.Simulation(
            # ------------------------------------------------------------
            # Domain: 200 x 100 cells, cell size 0.40 di in each direction.
            # Physical size: Lx = 80 di, Ly = 40 di.
            # The aspect ratio Lx/Ly = 2 is wide enough to host multiple
            # reconnection sites along x while keeping the run affordable.
            # ------------------------------------------------------------
            cells=(200, 100),
            dl=(0.40, 0.40),
            time_step=time_step,
            final_time=final_time,
            # ------------------------------------------------------------
            # Resistivity and hyper-resistivity control reconnection onset.
            # resistivity=0.001  provides the seed dissipation that breaks
            # ideal MHD; hyper_resistivity=0.002 damps grid-scale noise
            # without broadening the current sheet noticeably.
            # ------------------------------------------------------------
            hyper_resistivity=0.002,
            resistivity=0.001,
            # ------------------------------------------------------------
            # AMR: tagging-based refinement tracks the evolving current
            # sheets automatically. max_nbr_levels=1 keeps the run at two
            # grid levels (coarse + one refined level).
            # ------------------------------------------------------------
            refinement="tagging",
            max_nbr_levels=1,
            nesting_buffer=1,
            diag_options={
                "format": "phareh5",
                "options": {"dir": diag_dir, "mode": "overwrite"},
            },
            strict=True,
        )

        # Retrieve physical domain extents after the Simulation object is
        # built (they depend on cells * dl).
        Lx = sim.simulation_domain()[0]   # 80 di
        Ly = sim.simulation_domain()[1]   # 40 di

        # ------------------------------------------------------------
        # Double Harris density: two current sheets at y = 0.3 Ly and
        # y = 0.7 Ly, plus a uniform background density of 0.4.
        # ------------------------------------------------------------
        def density(x, y):
            return (
                0.4
                + 1.0 / np.cosh((y - Ly * 0.3) / L) ** 2
                + 1.0 / np.cosh((y - Ly * 0.7) / L) ** 2
            )

        # Smooth step function — used to build the Bx tanh profile below.
        def S(y, y0, l):
            return 0.5 * (1.0 + np.tanh((y - y0) / l))

        # ------------------------------------------------------------
        # Equilibrium Bx: reverses sign across each current sheet.
        # The S() terms produce +1 far from the sheet and -1 inside it.
        # A small Gaussian perturbation (amplitude dB = 0.1, width sigma)
        # seeds the tearing instability at both neutral lines.
        # ------------------------------------------------------------
        def bx(x, y):
            sigma = 1.0
            dB    = 0.1

            x0 = x - 0.5 * Lx
            y1 = y - 0.3 * Ly
            y2 = y - 0.7 * Ly

            dBx1 = -2 * dB * y1 * np.exp(-(x0**2 + y1**2) / sigma**2)
            dBx2 =  2 * dB * y2 * np.exp(-(x0**2 + y2**2) / sigma**2)

            return -1.0 + 2.0 * (S(y, Ly * 0.3, L) - S(y, Ly * 0.7, L)) + dBx1 + dBx2

        # ------------------------------------------------------------
        # By: zero in the equilibrium; the Gaussian perturbation adds a
        # small By component to seed the X-point geometry.
        # ------------------------------------------------------------
        def by(x, y):
            sigma = 1.0
            dB    = 0.1

            x0 = x - 0.5 * Lx
            y1 = y - 0.3 * Ly
            y2 = y - 0.7 * Ly

            dBy1 =  2 * dB * x0 * np.exp(-(x0**2 + y1**2) / sigma**2)
            dBy2 = -2 * dB * x0 * np.exp(-(x0**2 + y2**2) / sigma**2)

            return dBy1 + dBy2

        def bz(x, y):
            return 0.0

        # Total magnetic pressure — needed to compute the temperature.
        def b2(x, y):
            return bx(x, y) ** 2 + by(x, y) ** 2 + bz(x, y) ** 2

        # ------------------------------------------------------------
        # Temperature: derived from total pressure balance.
        # K = 0.7 is the total pressure constant.
        # The assertion checks that T > 0 everywhere (mandatory).
        # ------------------------------------------------------------
        def T(x, y):
            K    = 0.7
            temp = (K - 0.5 * b2(x, y)) / density(x, y)
            assert np.all(temp > 0)
            return temp

        # Bulk velocity: the Harris equilibrium is at rest.
        def vx(x, y):
            return 0.0

        def vy(x, y):
            return 0.0

        def vz(x, y):
            return 0.0

        # Thermal speed: isotropic Maxwellian with local temperature T(x,y).
        def vthx(x, y):
            return np.sqrt(T(x, y))

        def vthy(x, y):
            return np.sqrt(T(x, y))

        def vthz(x, y):
            return np.sqrt(T(x, y))

        # ------------------------------------------------------------
        # Particle model: single proton population, 100 macro-particles
        # per cell.  The seed fixes the random draw so runs are
        # reproducible.
        # ------------------------------------------------------------
        ph.MaxwellianFluidModel(
            bx=bx,
            by=by,
            bz=bz,
            protons={
                "charge": 1,
                "density": density,
                "vbulkx": vx,
                "vbulky": vy,
                "vbulkz": vz,
                "vthx": vthx,
                "vthy": vthy,
                "vthz": vthz,
                "nbr_part_per_cell": 100,
                "init": {"seed": 12334},
            },
        )

        # Cold isothermal electrons — removes the electron pressure
        # gradient from Ohm's law, isolating ion-scale reconnection physics.
        ph.ElectronModel(closure="isothermal", Te=0.0)

        # ------------------------------------------------------------
        # Diagnostics: electromagnetic fields and ion fluid moments.
        # ------------------------------------------------------------
        for quantity in ["E", "B"]:
            ph.ElectromagDiagnostics(quantity=quantity, write_timestamps=timestamps)

        for quantity in ["mass_density", "bulkVelocity"]:
            ph.FluidDiagnostics(quantity=quantity, write_timestamps=timestamps)

        for quantity in ["density", "pressure_tensor"]:
            ph.FluidDiagnostics(
                quantity=quantity,
                write_timestamps=timestamps,
                population_name="protons",
            )

        ph.InfoDiagnostics(quantity="particle_count")

        # Dynamic load balancing: redistribute patches across MPI ranks
        # based on particle count (nppc mode) when imbalance exceeds 5 %.
        ph.LoadBalancer(active=True, auto=True, mode="nppc", tol=0.05)

        return sim


    if __name__ == "__main__":
        startMPI()
        Simulator(config()).run()

For a complete description of every configuration block, see
:doc:`../usage/simulation_inputs`.


Key parameter choices
~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Parameter
     - Value
     - Rationale
   * - ``cells``
     - ``(200, 100)``
     - Aspect ratio 2:1 allows multiple X-points to form along x.
   * - ``dl``
     - ``(0.40, 0.40)``
     - Resolves the current sheet half-width *L* = 0.5 di with ~1.25 cells.
   * - ``L``
     - 0.5
     - Current sheet half-width; thin enough to be unstable, thick enough
       to be resolved.
   * - ``n_bg``
     - 0.4
     - Background density; keeps particle statistics acceptable away from
       the sheets.
   * - ``K``
     - 0.7
     - Total pressure constant; set so that *T* > 0 everywhere.
   * - ``dB``
     - 0.1
     - Perturbation amplitude; large enough to seed reconnection rapidly
       without pre-determining the X-point location.
   * - ``nbr_part_per_cell``
     - 100
     - Balances statistical noise against memory cost for a 2D run.
   * - ``resistivity``
     - 0.001
     - Provides the irreversible dissipation that triggers reconnection.
   * - ``hyper_resistivity``
     - 0.002
     - Damps grid-scale noise without broadening the macroscopic current
       sheet.


Running the simulation
-----------------------

2D simulations benefit significantly from MPI parallelism. Run with four
ranks for a workstation:

.. code-block:: bash

    mpirun -np 4 python3 -Ou harris_2d.py

The ``-O`` flag disables Python assertions for speed (remove it during
development to catch configuration errors early). The ``-u`` flag forces
unbuffered stdout so progress messages appear immediately.

On a cluster, submit a job script that passes the same command to your
scheduler. The ``LoadBalancer`` block rebalances work between ranks
automatically as the reconnection region concentrates particles.

Diagnostics accumulate in ``phare_outputs/harris_2d/`` as HDF5 files
(``EM_B.h5``, ``EM_E.h5``, ``ions_mass_density.h5``, etc.).


Analyzing results
-----------------

Load the output with ``pharesee`` and visualize the out-of-plane current
density :math:`J_z`, the clearest signature of active reconnection:

.. code-block:: python

    import numpy as np
    import matplotlib.pyplot as plt
    from pyphare.pharesee.run import Run
    from pyphare.pharesee.hierarchy.fromh5 import get_times_from_h5

    diag_dir = "phare_outputs/harris_2d"
    run = Run(diag_dir)

    times = get_times_from_h5(f"{diag_dir}/EM_B.h5")

    # Pick three snapshots: beginning, mid-run, and end.
    selected = [times[0], times[len(times) // 2], times[-1]]

    fig, axes = plt.subplots(1, 3, figsize=(15, 4), constrained_layout=True)

    for ax, t in zip(axes, selected):
        J = run.GetJ(t, merged=True)
        # GetJ returns (interpolant, [x_coords, y_coords]) for each component.
        xi, yi = J["Jz"][1]
        Jz = J["Jz"][0](xi, yi)

        pcm = ax.pcolormesh(xi, yi, Jz.T, cmap="RdBu_r", vmin=-2, vmax=2)
        ax.set_title(f"$J_z$,  t = {t:.1f}")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        fig.colorbar(pcm, ax=ax, label=r"$J_z$")

    plt.show()

Reconnection onset is visible as a strong positive or negative :math:`J_z`
filament at the X-point, flanked by two magnetic islands (O-points) with
opposite :math:`J_z` sign.


Reconnection rate analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~

The reconnection electric field (reconnection rate) is computed from the
time derivative of the magnetic flux :math:`A_z` at the X-point. ``pharesee``
provides dedicated helpers:

.. code-block:: python

    import numpy as np
    import matplotlib.pyplot as plt
    from pyphare.pharesee.run import Run
    from pyphare.pharesee.hierarchy.fromh5 import get_times_from_h5

    diag_dir = "phare_outputs/harris_2d"
    run = Run(diag_dir)

    times = get_times_from_h5(f"{diag_dir}/EM_B.h5")

    # Magnetic flux function Az at the last snapshot.
    Az, (xn, yn) = run.GetMagneticFlux(times[-1])

    # Reconnection rate time series: centered finite differences of Az at
    # the X-point location, returned together with the X-point trajectory.
    times_centered, rates, flux_at_xpoint, xpoint_traj = run.GetReconnectionRate(times)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4), constrained_layout=True)

    axes[0].pcolormesh(xn, yn, Az.T, cmap="RdBu_r")
    axes[0].plot(*xpoint_traj[-1], "k+", markersize=12, label="X-point")
    axes[0].set_title(f"$A_z$ at t = {times[-1]:.1f}")
    axes[0].set_xlabel("x")
    axes[0].set_ylabel("y")
    axes[0].legend()

    axes[1].plot(times_centered, rates)
    axes[1].set_xlabel("time")
    axes[1].set_ylabel("reconnection rate")
    axes[1].set_title("Reconnection electric field at X-point")

    plt.show()

``GetMagneticFlux(time)`` integrates :math:`B_y` along x to recover
:math:`A_z(x, y)`. ``GetReconnectionRate(times)`` locates the X-point at
each time (minimum of :math:`|A_z|` on the neutral line), evaluates
:math:`A_z` there, and differentiates with respect to time.

For a complete reference on these utilities, see
:doc:`../pharesee/reconnection_analysis`.


Expected signatures
~~~~~~~~~~~~~~~~~~~

A healthy reconnection run shows:

- :math:`J_z` filaments forming at both neutral lines (at :math:`y \approx 12`
  and :math:`y \approx 28` for the default domain).
- :math:`A_z` contours developing an X-point topology; closed contours
  (magnetic islands) growing symmetrically on each side.
- The reconnection rate rising from zero, peaking near the nonlinear
  saturation of the tearing mode, then declining as the islands merge.

.. note::

   With ``final_time = 0.005`` (the value in the regression test) only a
   single time step runs — enough to verify the code executes without error
   but not enough to observe reconnection. Set ``final_time = 50.0`` and
   choose meaningful ``timestamps`` to reproduce the physics described above.


Next steps
----------

- Increase ``max_nbr_levels`` to 2 or 3 to follow the thin current layer
  at the X-point with adaptively refined cells while keeping coarse resolution
  elsewhere.
- Add a second ion population (e.g. a cold background beam) to study
  reconnection in multi-fluid plasmas.
- Enable ``restart_options`` to checkpoint long runs and resume from any
  saved time.
- Compare the reconnection rate with resistive MHD predictions to quantify
  the role of kinetic effects resolved by the PIC ions.

See also:

- :doc:`../usage/simulation_inputs` — full reference for all configuration blocks
- :doc:`../theory/amr` — how PHARE's adaptive mesh refinement tracks evolving
  plasma structures
- :doc:`../pharesee/reconnection_analysis` — detailed API for flux and
  reconnection rate diagnostics
- :doc:`../pharesee/get_data` — general guide to reading simulation output
