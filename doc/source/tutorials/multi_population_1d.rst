
Tutorial: Multi-Population Ion-Ion Beam
========================================

This tutorial demonstrates PHARE's multi-population capability by simulating
two counter-streaming proton beams in 1D. It introduces how to configure
multiple ion populations with distinct bulk velocities and how to diagnose
each population independently. By the end you will have measured the
electromagnetic growth rate of the beam instability and produced phase-space
plots showing how the beam structure evolves.

.. contents:: Contents
   :local:
   :depth: 2


Introduction
------------

Counter-streaming ion beams are ubiquitous in space plasmas: the solar wind
carries a secondary proton beam that propagates along the interplanetary
magnetic field at roughly the Alfven speed, and magnetosheath jets contain
fast ion populations that interpenetrate the ambient plasma. When two ion
populations drift relative to each other along the background magnetic field,
the free energy in the relative drift drives electromagnetic instabilities.

This scenario also serves as a standard quantitative benchmark: for a
specific set of parameters the linear theory predicts a right-hand
circularly polarized wave growing at a well-defined rate
:math:`\gamma \approx 0.09 \,\Omega_{ci}`, and the simulation can be
verified against this prediction.


The physics
-----------

Two proton populations share the same 1D periodic box:

- **Main population** — at rest (:math:`V_{b,main} = 0`), density
  :math:`n_{main} = 1` (normalized).
- **Beam population** — drifting at :math:`V_{b,beam} = 5\,V_A` along
  *x*, density :math:`n_{beam} = 0.01 \, n_{main}`.

A uniform background magnetic field :math:`B_0 = 1` is oriented along *x*.
Both populations have an isotropic thermal speed
:math:`v_{th} = \sqrt{0.1} \approx 0.316\,V_A`.

The relative drift excites the right-hand resonant beam instability
(sometimes called the ion/ion right-hand resonant instability). Linear
theory for these parameters predicts the most unstable mode at
:math:`k d_i \approx 0.19`, corresponding to a wavelength
:math:`\lambda \approx 33\,d_i` (ion inertial lengths). The domain length
is chosen to be exactly one wavelength so that the Fourier mode :math:`m = 1`
captures the most unstable perturbation.

The instability grows exponentially during the linear phase, saturates, and
then enters a nonlinear phase where the beam loses coherence and the
populations mix in phase space.

For theoretical background on the hybrid PIC model, see
:doc:`../theory/hybridpic`.


Setting up the simulation
--------------------------

Create a file ``ion_ion_beam1d.py`` with the content below. The key new
element compared to a single-population run is the ``MaxwellianFluidModel``
call: population names are arbitrary Python keyword arguments, so adding a
second population is as simple as adding a second keyword.

.. code-block:: python

    #!/usr/bin/env python3

    import numpy as np
    import pyphare.pharein as ph
    from pyphare.simulator.simulator import Simulator

    ph.NO_GUI()


    def config():
        # ------------------------------------------------------------------
        # Domain: 165 cells x dl = 0.2 => physical length L = 33 d_i.
        # This fits exactly one wavelength of the most unstable mode at
        # k d_i = 0.19  (lambda = 2*pi/k ~ 33).
        # ------------------------------------------------------------------
        sim = ph.Simulation(
            smallest_patch_size=20,
            largest_patch_size=60,
            final_time=100,
            time_step=0.001,
            boundary_types="periodic",
            cells=165,
            dl=0.2,                       # cell size in ion inertial lengths
            hyper_resistivity=0.01,       # damps sub-grid noise
            refinement="tagging",         # AMR: refine where gradients are large
            max_nbr_levels=3,
            diag_options={
                "format": "phareh5",
                "options": {"dir": "ion_ion_beam1d", "mode": "overwrite"},
            },
        )

        # ------------------------------------------------------------------
        # Background field: B0 = 1 along x (sets VA = 1 in normalized units).
        # By = Bz = 0 initially; perturbations grow from numerical noise.
        # ------------------------------------------------------------------
        def bx(x):
            return 1.0

        def by(x):
            return 0.0

        def bz(x):
            return 0.0

        # ------------------------------------------------------------------
        # Main population: at rest, density = 1.
        # ------------------------------------------------------------------
        def density_main(x):
            return 1.0

        # ------------------------------------------------------------------
        # Beam population: drifting at vB = 5 VA along x, density = 0.01.
        # A low density ensures the beam is a minor component (realistic for
        # the solar wind secondary beam).
        # ------------------------------------------------------------------
        def density_beam(x):
            return 0.01

        def v_beam(x):
            return 5.0

        def v_zero(x):
            return 0.0

        # Thermal speed: vth = sqrt(0.1) ~ 0.316 VA for both populations.
        def vth(x):
            return np.sqrt(0.1)

        # ------------------------------------------------------------------
        # MaxwellianFluidModel: each keyword argument after bx/by/bz defines
        # one ion population. The name (here "main" and "beam") is used to
        # label diagnostics files and is entirely user-chosen.
        # ------------------------------------------------------------------
        ph.MaxwellianFluidModel(
            bx=bx,
            by=by,
            bz=bz,
            main={
                "charge": 1,
                "density": density_main,
                "vbulkx": v_zero,
                "vbulky": v_zero,
                "vbulkz": v_zero,
                "vthx": vth,
                "vthy": vth,
                "vthz": vth,
            },
            beam={
                "charge": 1,
                "density": density_beam,
                "vbulkx": v_beam,   # drift along B0
                "vbulky": v_zero,
                "vbulkz": v_zero,
                "vthx": vth,
                "vthy": vth,
                "vthz": vth,
            },
        )

        # ------------------------------------------------------------------
        # Electron model: cold isothermal electrons (Te = 0).
        # ------------------------------------------------------------------
        ph.ElectronModel(closure="isothermal", Te=0.0)

        # ------------------------------------------------------------------
        # Diagnostics: write snapshots every 0.1 ion cyclotron periods.
        # Per-population particle diagnostics allow tracking each beam
        # independently in phase space.
        # ------------------------------------------------------------------
        timestamps = np.arange(0, sim.final_time, 0.1)

        for quantity in ["B", "E"]:
            ph.ElectromagDiagnostics(quantity=quantity, write_timestamps=timestamps)

        for pop_name in ["main", "beam"]:
            ph.ParticleDiagnostics(
                quantity="domain",
                population_name=pop_name,
                write_timestamps=timestamps,
            )

        return sim


    if __name__ == "__main__":
        Simulator(config()).run()

.. note::

   Population names passed to ``MaxwellianFluidModel`` are arbitrary Python
   keyword arguments. You can add a third population (e.g. ``alphas={...}``)
   without changing anything else. The name is propagated to the HDF5 output
   files so each population's moments and particles are stored separately.

For the full reference of every accepted parameter, see
:doc:`../usage/simulation_inputs`.


Running the simulation
-----------------------

From the directory containing ``ion_ion_beam1d.py``, run:

.. code-block:: bash

    python3 -Ou ion_ion_beam1d.py

The ``-O`` flag disables Python assertions (faster) and ``-u`` forces
unbuffered output. On a modern laptop this run completes in roughly 15–30
minutes because the domain is small (165 cells, :math:`\Delta t = 0.001`) and
the final time is 100 ion cyclotron periods.

To accelerate the run with MPI:

.. code-block:: bash

    mpirun -n 4 python3 -Ou ion_ion_beam1d.py

Output files are written to the subdirectory ``ion_ion_beam1d/``:

- ``EM_B.h5`` — magnetic field snapshots
- ``EM_E.h5`` — electric field snapshots
- ``ions_pop_main_domain.h5`` — main population particles
- ``ions_pop_beam_domain.h5`` — beam population particles


Analyzing results
-----------------

Density profiles
~~~~~~~~~~~~~~~~

Compare the density of each population as a function of position at several
times. Initially both are uniform; once the instability saturates the beam
density develops spatial structure.

.. code-block:: python

    import numpy as np
    import matplotlib.pyplot as plt
    from pyphare.pharesee.run import Run
    from pyphare.pharesee.hierarchy.fromh5 import get_times_from_h5

    run_path = "ion_ion_beam1d"
    r = Run(run_path)

    times = get_times_from_h5(f"{run_path}/EM_B.h5")
    t_early = times[0]
    t_peak  = times[len(times) // 2]
    t_late  = times[-1]

    fig, axes = plt.subplots(1, 3, figsize=(12, 4), sharey=True)

    for ax, t in zip(axes, [t_early, t_peak, t_late]):
        for pop in ["main", "beam"]:
            ions = r.GetParticles(t, [pop])
            # GetParticles returns a MultiLevelParticleData object;
            # dist_plot projects particles onto a 1D histogram
            ions.dist_plot(
                axis=("x",),
                ax=ax,
                label=pop,
            )
        ax.set_title(f"t = {t:.1f}")
        ax.set_xlabel("x")

    axes[0].set_ylabel("density (arb.)")
    axes[0].legend()
    fig.tight_layout()
    plt.show()

Bulk velocity of each population
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Track the bulk velocity profiles to observe how the beam decelerates and the
main population accelerates as momentum is exchanged via the wave.

.. code-block:: python

    import numpy as np
    import matplotlib.pyplot as plt
    from pyphare.pharesee.run import Run
    from pyphare.pharesee.hierarchy.fromh5 import get_times_from_h5

    run_path = "ion_ion_beam1d"
    r = Run(run_path)
    times = get_times_from_h5(f"{run_path}/EM_B.h5")

    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)

    for ax, pop in zip(axes, ["main", "beam"]):
        for t in [times[0], times[len(times) // 2], times[-1]]:
            ions = r.GetParticles(t, [pop])
            ions.dist_plot(
                axis=("x",),
                moment="vx",
                ax=ax,
                label=f"t = {t:.0f}",
            )
        ax.set_title(pop)
        ax.set_xlabel("x")
        ax.legend()

    axes[0].set_ylabel(r"$V_{bulk,x}$")
    fig.tight_layout()
    plt.show()

Phase-space plot (x vs Vx)
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The most revealing diagnostic for beam instabilities is a joint
(x, Vx) phase-space plot. Initially the two populations appear as two
horizontal bands (main near :math:`V_x = 0`, beam near :math:`V_x = 5`).
As the instability grows, the beam develops sinusoidal oscillations in Vx
(phase bunching), and at saturation the particles form vortices (holes) in
phase space — a signature of trapping.

.. code-block:: python

    import numpy as np
    import matplotlib.pyplot as plt
    from pyphare.pharesee.run import Run
    from pyphare.pharesee.hierarchy.fromh5 import get_times_from_h5

    run_path = "ion_ion_beam1d"
    r = Run(run_path)
    times = get_times_from_h5(f"{run_path}/EM_B.h5")

    # Pick three representative times: initial, peak instability, nonlinear
    t_early = times[0]
    t_peak  = times[len(times) // 2]
    t_late  = times[-1]

    fig, axes = plt.subplots(1, 3, figsize=(13, 5), constrained_layout=True)
    vmin, vmax = -1.3, 5.6

    for ax, t in zip(axes, [t_early, t_peak, t_late]):
        ions = r.GetParticles(t, ["main", "beam"])
        ions.dist_plot(
            axis=("x", "Vx"),
            ax=ax,
            norm=0.4,
            finest=True,
            gaussian_filter_sigma=(1, 1),
            vmin=vmin,
            vmax=vmax,
            dv=0.05,
            title=f"t = {t:.1f}",
            xlabel="x",
            ylabel=r"$V_x$",
        )
        ax.set_xlim((0, 33))

    fig.savefig("phase_space.png", dpi=150)
    plt.show()

.. note::

   ``r.GetParticles(t, ["main", "beam"])`` combines both populations into a
   single ``MultiLevelParticleData`` object for a joint distribution. To plot
   them separately, call ``r.GetParticles(t, ["main"])`` and
   ``r.GetParticles(t, ["beam"])`` individually.

Measuring the growth rate
~~~~~~~~~~~~~~~~~~~~~~~~~

The benchmark quantity is the growth rate :math:`\gamma` of the right-hand
circular polarization mode at :math:`m = 1`. It is extracted from the
time series of the first Fourier mode of the complex field
:math:`B_y - i B_z`:

.. code-block:: python

    import os
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.optimize import curve_fit
    from scipy.signal import find_peaks
    from pyphare.pharesee.run import Run
    from pyphare.pharesee.hierarchy.fromh5 import get_times_from_h5

    run_path = "ion_ion_beam1d"
    r = Run(run_path)
    times = get_times_from_h5(os.path.join(run_path, "EM_B.h5"))

    first_mode = []
    for t in times:
        B = r.GetB(t, merged=True, interp="linear")
        by_interp, xyz = B["By"]
        bz_interp, _   = B["Bz"]
        x = xyz[0][:-1]                         # drop last (= first) point for periodicity
        by = by_interp(x)
        bz = bz_interp(x)
        mode1 = np.abs(np.fft.fft(by - 1j * bz)[1])
        first_mode.append(mode1)

    first_mode = np.array(first_mode)
    dt = times[1] - times[0]

    # Locate the peak of the instability
    time_offset = 10.0
    imax = find_peaks(first_mode, width=int(time_offset / dt))[0][0]

    # Fit an exponential to the linear phase
    def exp_growth(t, a, gamma):
        return a * np.exp(gamma * t)

    popt, _ = curve_fit(
        exp_growth,
        times[: imax - int(time_offset / dt)],
        first_mode[: imax - int(time_offset / dt)],
        p0=[0.08, 0.09],
    )
    gamma = popt[1]
    print(f"Measured growth rate: gamma = {gamma:.4f}  (expected ~0.09)")

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.semilogy(times, first_mode, "k", label=r"$|\hat{B}_y(m=1) - i\hat{B}_z(m=1)|$")
    ax.semilogy(
        times[:imax],
        exp_growth(times[:imax], *popt),
        "r--",
        label=rf"$B_0 \exp(\gamma t)$, $\gamma = {gamma:.4f}$",
    )
    ax.set_xlabel(r"$t \; [\Omega_{ci}^{-1}]$")
    ax.set_ylabel("mode amplitude")
    ax.legend()
    fig.tight_layout()
    fig.savefig("growth_rate.png", dpi=150)
    plt.show()

A correctly configured run yields :math:`\gamma` within about 2% of the
theoretical value :math:`\gamma_{th} \approx 0.09\,\Omega_{ci}`.


Expected behavior
-----------------

- **t = 0** — both populations are spatially uniform; the magnetic field is
  purely along *x*. The by and bz components start at zero and grow from
  numerical noise.
- **Linear phase (t ≈ 0–60)** — the right-hand circular mode at
  :math:`m = 1` grows exponentially at rate :math:`\gamma \approx 0.09`.
  The beam particles develop a sinusoidal modulation in Vx.
- **Saturation (t ≈ 60–80)** — wave amplitude peaks; particles form
  closed loops (holes) in (x, Vx) phase space — the signature of particle
  trapping.
- **Nonlinear phase (t > 80)** — the beam loses its coherent drift structure;
  energy is redistributed among the wave modes and both populations.

.. note::

   AMR is enabled via ``refinement="tagging"`` and ``max_nbr_levels=3``.
   PHARE will automatically refine the mesh in regions where electromagnetic
   gradients are large, concentrating resolution where the instability is
   most active. You can inspect the refined patches by examining the patch
   hierarchy in the HDF5 output.


Next steps
----------

- Add fluid diagnostics (``FluidDiagnostics`` with ``quantity="density"``
  and ``quantity="bulkVelocity"``) for each population to access
  macroscopic moments without loading all particles.
- Increase ``max_nbr_levels`` to 4 or 5 and compare the growth rate with
  the uniform-grid result to assess the AMR overhead.
- Replace the isothermal electron closure with ``"adiabatic"`` and a finite
  electron temperature to study how :math:`T_e` modifies the instability
  threshold.
- Run the simulation in 2D by switching to a 2D ``cells`` tuple and
  observing the formation of obliquely propagating modes.

See also:

- :doc:`../usage/simulation_inputs` — full reference for all configuration blocks
- :doc:`../pharesee/get_data` — reading and manipulating simulation output
- :doc:`../pharesee/plotting_distributions` — phase-space visualization utilities
- :doc:`../theory/hybridpic` — theoretical background of the hybrid PIC model
