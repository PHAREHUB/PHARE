
Tutorial: 1D Alfven Wave
========================

This tutorial walks through your first PHARE simulation: a circularly polarized
Alfven wave propagating in a uniform magnetized plasma. It introduces the core
building blocks of any PHARE input script and shows how to verify the result
by measuring the wave phase speed from the output data.

.. contents:: Contents
   :local:
   :depth: 2


The physics
-----------

An Alfven wave is a transverse magnetohydrodynamic wave that propagates along
the background magnetic field **B**\ :sub:`0`. In a uniform plasma of density
*n*\ :sub:`0` and ion mass *m*, the Alfven speed is:

.. math::

   V_A = \frac{B_0}{\sqrt{\mu_0 \, n_0 \, m}}

In PHARE's normalized units (where :math:`\mu_0 = 1`, *m* = 1, and we choose
:math:`B_0 = n_0 = 1`), the Alfven speed is simply :math:`V_A = 1`.

A circularly polarized Alfven wave consists of sinusoidal perturbations in the
transverse components of both the magnetic field and the ion bulk velocity,
phase-locked so that:

.. math::

   \delta B_y = \delta B \cos(kx), \quad \delta B_z = \delta B \sin(kx)

   \delta V_y = \delta V \cos(kx), \quad \delta V_z = \delta V \sin(kx)

For a forward-propagating wave, the velocity and magnetic perturbations satisfy
:math:`\delta V / V_A = \delta B / B_0`. This test is a standard validation
benchmark for hybrid PIC codes: after a simulation time *T*, the wave pattern
should have shifted by exactly :math:`V_A \times T` in the +x direction.

For background reading on the hybrid PIC model, see :doc:`../theory/hybridpic`.


Setting up the simulation
--------------------------

Create a file ``alfven_wave_1d.py`` with the content below.
Each parameter choice is explained in the comments.

.. code-block:: python

    #!/usr/bin/env python3

    import numpy as np
    import pyphare.pharein as ph
    from pyphare.simulator.simulator import Simulator

    ph.NO_GUI()


    def config():
        # ------------------------------------------------------------
        # Domain: 1000 cells of size dl=1 => physical length L = 1000.
        # This fits exactly one wavelength of a k = 2*pi/L perturbation.
        # ------------------------------------------------------------
        final_time = 1000          # VA = 1, so the wave travels L = 1000
        time_step_nbr = 100000     # dt = final_time / time_step_nbr = 0.01

        sim = ph.Simulation(
            smallest_patch_size=50,
            largest_patch_size=50,
            time_step_nbr=time_step_nbr,
            final_time=final_time,
            boundary_types="periodic",  # wave re-enters from the left
            cells=1000,                 # number of grid cells (1D)
            dl=1,                       # cell size in ion inertial lengths
            # hyper-resistivity damps grid-scale noise without affecting
            # the long-wavelength wave we are studying
            hyper_resistivity=0.001,
            diag_options={
                "format": "phareh5",
                "options": {"dir": ".", "mode": "overwrite"},
            },
        )

        # L is the physical domain length (= cells * dl = 1000)
        L = sim.simulation_domain()[0]

        # ------------------------------------------------------------
        # Background field: B0 = 1 along x (sets VA = 1 in norm. units)
        # ------------------------------------------------------------
        def bx(x):
            return 1.0

        # ------------------------------------------------------------
        # Perturbation: amplitude 0.01, one full wavelength across L.
        # By and Bz form a right-hand circularly polarized wave.
        # ------------------------------------------------------------
        def by(x):
            return 0.01 * np.cos(2 * np.pi * x / L)

        def bz(x):
            return 0.01 * np.sin(2 * np.pi * x / L)

        # ------------------------------------------------------------
        # Ion density: uniform (n0 = 1)
        # ------------------------------------------------------------
        def density(x):
            return 1.0

        # ------------------------------------------------------------
        # Ion bulk velocity: no mean flow (vx = 0).
        # vy and vz match the magnetic perturbation amplitudes so that
        # delta_V / VA = delta_B / B0  (Alfven wave dispersion relation).
        # With VA = B0 = 1 and delta_B = 0.01 this gives delta_V = 0.01.
        # ------------------------------------------------------------
        def vx(x):
            return 0.0

        def vy(x):
            return 0.01 * np.cos(2 * np.pi * x / L)

        def vz(x):
            return 0.01 * np.sin(2 * np.pi * x / L)

        # ------------------------------------------------------------
        # Thermal velocity: cold plasma (vth << VA and << delta_V).
        # A small but nonzero value avoids singular particle distributions.
        # ------------------------------------------------------------
        def vthx(x):
            return 0.01

        def vthy(x):
            return 0.01

        def vthz(x):
            return 0.01

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
            },
        )

        # ------------------------------------------------------------
        # Electron model: cold isothermal electrons (Te = 0).
        # Setting Te = 0 removes the electron pressure gradient from
        # the generalized Ohm's law, keeping the setup as close as
        # possible to ideal MHD for this wave test.
        # ------------------------------------------------------------
        ph.ElectronModel(closure="isothermal", Te=0.0)

        # ------------------------------------------------------------
        # Diagnostics: write 11 snapshots evenly spaced in time
        # (t = 0, 100, 200, ..., 1000) for both E and B fields and
        # for the ion fluid moments.
        # ------------------------------------------------------------
        timestamps = np.arange(0, sim.final_time + sim.time_step, sim.final_time / 10)

        for quantity in ["E", "B"]:
            ph.ElectromagDiagnostics(quantity=quantity, write_timestamps=timestamps)

        for quantity in ["charge_density", "bulkVelocity"]:
            ph.FluidDiagnostics(quantity=quantity, write_timestamps=timestamps)

        return sim


    if __name__ == "__main__":
        Simulator(config()).run()

For details on every parameter accepted by each block, see
:doc:`../usage/simulation_inputs`.


Running the simulation
-----------------------

From the directory containing ``alfven_wave_1d.py``, run:

.. code-block:: bash

    python3 -Ou alfven_wave_1d.py

The ``-O`` flag disables assertion checks (faster execution) and ``-u``
forces unbuffered output so progress messages appear immediately. On a
laptop the run takes a few minutes. The diagnostics are written to the
current directory (``dir: "."``) as HDF5 files such as ``EM_B.h5``.

To use multiple MPI ranks (e.g. on a workstation or cluster):

.. code-block:: bash

    mpirun -n 4 python3 -Ou alfven_wave_1d.py


Analyzing results
-----------------

After the simulation completes, use ``pharesee`` to read the output and
plot :math:`B_y` at several times.

.. code-block:: python

    import numpy as np
    import matplotlib.pyplot as plt
    from pyphare.pharesee.run import Run
    from pyphare.pharesee.hierarchy.fromh5 import get_times_from_h5

    run = Run(".")   # path to the directory containing EM_B.h5

    times = get_times_from_h5("EM_B.h5")

    fig, ax = plt.subplots(figsize=(9, 4))

    for t in [times[0], times[len(times) // 2], times[-1]]:
        B = run.GetB(t, merged=True)
        x = B["By"][1][0]       # coordinate array
        by = B["By"][0](x)      # By values interpolated onto x
        ax.plot(x, by, label=f"t = {t:.0f}")

    ax.set_xlabel("x")
    ax.set_ylabel(r"$B_y$")
    ax.legend()
    fig.tight_layout()
    plt.show()

``run.GetB(t, merged=True)`` returns a dictionary whose keys are
``"Bx"``, ``"By"``, ``"Bz"``. Each value is a two-element tuple
``(interpolant, [coordinate_array])``. Calling ``interpolant(x)``
evaluates the field on any coordinate array you supply.

For a complete guide to reading simulation output, see
:doc:`../pharesee/get_data`.


Verifying the Alfven speed
--------------------------

The wave should propagate at exactly :math:`V_A = 1`. One way to verify
this is to track the phase of the :math:`B_y` sinusoid over time using
a curve fit, then compute :math:`d\phi/dt \;/\; k`.

.. code-block:: python

    import numpy as np
    from scipy.signal import medfilt
    from scipy.optimize import curve_fit
    from pyphare.pharesee.run import Run
    from pyphare.pharesee.hierarchy.fromh5 import get_times_from_h5


    def wave(x, a0, k, phi):
        return a0 * np.cos(k * x + phi)


    run_path = "."
    time = get_times_from_h5("EM_B.h5")
    r = Run(run_path)

    L = 1000.0
    phase = np.zeros_like(time)
    wave_vec = np.zeros_like(time)

    for it, t in enumerate(time):
        B = r.GetB(t, merged=True)
        x = B["By"][1][0]
        by = B["By"][0](x)
        a, k, phi = curve_fit(wave, x, by, p0=(0.01, 2 * np.pi / L, 0))[0]
        phase[it] = phi
        wave_vec[it] = k

    # Phase speed: d(phi)/dt / k. Median filter removes outliers.
    vphi = medfilt(np.gradient(phase, time) / wave_vec, kernel_size=7)
    print(f"Mean phase speed: {vphi.mean():.4f}  (expected 1.0000)")

A correctly configured run yields a phase speed within 5% of 1.0 — the
measured deviation reflects finite particle statistics and the small
hyper-resistive damping chosen for this test.


Expected output
---------------

At :math:`t = 0` the :math:`B_y` profile is a pure cosine across the domain.
By :math:`t = 500` (half an Alfven crossing time) the pattern should have
shifted by 500 grid cells, and by :math:`t = 1000` (one full crossing) it
returns to a shape close to the initial condition (with modest amplitude
reduction from hyper-resistivity).

.. note::

   The simulation uses periodic boundary conditions, so a phase shift
   larger than the domain length *L* simply wraps around. The curve-fit
   approach above tracks the cumulative phase and is not affected by this.


Next steps
----------

- Increase the perturbation amplitude to study nonlinear effects.
- Add a second ion population (e.g. alpha particles) to observe wave
  damping through resonant wave-particle interactions.
- Enable adaptive mesh refinement via ``refinement="tagging"`` and compare
  with the uniform-grid result.

See also:

- :doc:`../usage/simulation_inputs` — full reference for all configuration blocks
- :doc:`../pharesee/get_data` — reading and manipulating simulation output
- :doc:`../theory/hybridpic` — theoretical background of the hybrid PIC model
