
Tutorial: 1D MHD Shock Tube
============================

This tutorial walks through your first MHD simulation in PHARE: a one-dimensional
shock tube (Riemann problem). The shock tube is the fundamental validation test for
any MHD solver — it has a known analytical structure and exercises all the wave
families the equations support. It is a natural starting point before moving on to
multi-dimensional or adaptive-mesh MHD runs.

.. contents:: Contents
   :local:
   :depth: 2


The physics
-----------

A shock tube initializes the domain with two uniform states separated by a
discontinuity at the mid-point. When the simulation starts the discontinuity
breaks up into a set of waves that propagate away from each other. For ideal
gas dynamics (Sod problem) there are three waves: a rarefaction fan moving left,
a contact discontinuity, and a shock moving right.

Magnetohydrodynamics enriches this picture because the magnetic field introduces
two additional wave families. The seven MHD characteristics are:

- **fast magnetosonic** (rightward and leftward)
- **Alfven** (rightward and leftward)
- **slow magnetosonic** (rightward and leftward)
- **entropy / contact** wave

The relative ordering and amplitude of these waves depend on the ratio of the
magnetic pressure to the thermal pressure (plasma beta) and on the orientation of
**B** relative to the shock normal.

This test uses the Brio–Wu initial conditions, which are a standard MHD analogue
of the Sod problem. The left state is dense and highly pressurized; the right state
is tenuous and cold. A guide field :math:`B_x = 0.75` is constant everywhere.
The transverse field :math:`B_y` reverses sign across the interface, so the
configuration is a current sheet that is also a pressure discontinuity.

The initial left and right states are:

.. list-table::
   :header-rows: 1
   :widths: 20 20 20

   * - Quantity
     - Left state (:math:`x < L/2`)
     - Right state (:math:`x \geq L/2`)
   * - :math:`\rho`
     - 1.0
     - 0.125
   * - :math:`v_x, v_y, v_z`
     - 0, 0, 0
     - 0, 0, 0
   * - :math:`B_x`
     - 0.75
     - 0.75
   * - :math:`B_y`
     - 1.0
     - -1.0
   * - :math:`B_z`
     - 0.0
     - 0.0
   * - :math:`P`
     - 1.0
     - 0.1

All quantities are in PHARE's dimensionless units. The adiabatic index is
:math:`\gamma = 5/3`.

For theoretical background on the MHD model, see :doc:`../theory/mhd`.


Setting up the simulation
--------------------------

Create a file ``mhd_shock_1d.py`` with the content below. Each parameter choice is
explained in the comments.

.. code-block:: python

    #!/usr/bin/env python3

    import numpy as np
    import pyphare.pharein as ph
    from pyphare.simulator.simulator import Simulator, startMPI

    ph.NO_GUI()

    final_time  = 80      # run long enough for waves to separate clearly
    time_step   = 0.2     # explicit time step; keep dt small for stability
    timestamps  = [final_time]
    diag_dir    = "phare_outputs/shock"


    def config():
        # ----------------------------------------------------------------
        # Domain: 800 cells of size dl = 1  =>  physical length L = 800.
        # The discontinuity sits at x = 400 (L/2).
        # ----------------------------------------------------------------
        cells = (800,)
        dl    = (1.0,)

        sim = ph.Simulation(
            smallest_patch_size=15,
            time_step=time_step,
            final_time=final_time,
            cells=cells,
            dl=dl,
            # --- MHD-specific parameters ---
            # Select the MHD fluid model instead of the default hybrid-PIC model.
            model_options=["MHDModel"],
            # Reconstruction scheme for interface states.
            # "constant" is first-order (piecewise-constant); use "linear" with
            # a limiter for second-order accuracy.
            reconstruction="constant",
            limiter="",
            # Riemann solver used to compute numerical fluxes at cell interfaces.
            # "rusanov" is the simplest (Local Lax–Friedrichs); it is diffusive
            # but robust and appropriate for a first run.
            riemann="rusanov",
            # Time integration scheme.  "euler" is forward Euler (first order).
            mhd_timestepper="euler",
            # Adiabatic index gamma = 5/3 for a monatomic ideal gas.
            gamma=5.0 / 3.0,
            # Dissipation: set both resistivity and hyper-resistivity to zero
            # for an ideal MHD run.
            resistivity=0.0,
            eta=0.0,
            hyper_resistivity=0.0,
            nu=0.0,
            # AMR: single level (no refinement) for this introductory run.
            refinement="tagging",
            max_mhd_level=1,
            max_nbr_levels=1,
            diag_options={
                "format": "phareh5",
                "options": {"dir": diag_dir, "mode": "overwrite"},
            },
            strict=True,
        )

        # ----------------------------------------------------------------
        # Initial conditions: piecewise-constant left / right states.
        # np.where selects the left value where x < L/2, right value elsewhere.
        # ----------------------------------------------------------------
        L_half = cells[0] * dl[0] / 2   # = 400.0

        def density(x):
            return np.where(x < L_half, 1.0, 0.125)

        def vx(x): return 0.0
        def vy(x): return 0.0
        def vz(x): return 0.0

        # Guide field: constant Bx = 0.75 (normal to the discontinuity surface).
        def bx(x):
            return 0.75

        # Transverse field: reverses sign across the interface.
        def by(x):
            return np.where(x < L_half, 1.0, -1.0)

        def bz(x): return 0.0

        def p(x):
            return np.where(x < L_half, 1.0, 0.1)

        # Register the MHD model with PHARE.  All fluid quantities are
        # required: density, three velocity components, three B components,
        # and thermal pressure.
        ph.MHDModel(density=density, vx=vx, vy=vy, vz=vz,
                    bx=bx, by=by, bz=bz, p=p)

        # ----------------------------------------------------------------
        # Diagnostics: write selected MHD fluid quantities at final time.
        # ----------------------------------------------------------------
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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``model_options=["MHDModel"]``
    Switches PHARE from the default hybrid kinetic-ion / fluid-electron model
    to a fully fluid MHD model. Without this option PHARE expects particle
    populations defined with ``MaxwellianFluidModel``.

``reconstruction="constant"``
    First-order spatial reconstruction (piecewise-constant states at each cell
    face). Simple and robust, at the cost of significant numerical diffusion.
    For sharper wave structures, switch to ``"linear"`` and choose a limiter
    such as ``"minmod"`` or ``"vanleer"``.

``riemann="rusanov"``
    The Rusanov (Local Lax–Friedrichs) approximate Riemann solver. It uses the
    maximum signal speed in each cell to set a dissipative flux. More accurate
    solvers such as ``"hlld"`` are available and produce cleaner wave structures,
    but Rusanov is the safest choice when first testing a new configuration.

``mhd_timestepper="euler"``
    First-order explicit (forward Euler) time integration. Stable provided the
    CFL condition is satisfied. Higher-order options are available for production
    runs.

``gamma=5.0/3.0``
    Adiabatic index. :math:`5/3` is the standard value for a monatomic ideal gas
    and is used by the Brio–Wu reference solution.

``resistivity``, ``eta``, ``hyper_resistivity``, ``nu``
    All set to zero for ideal MHD. Non-zero values introduce Ohmic and viscous
    dissipation. They are exposed explicitly here to make the ideal-MHD assumption
    visible.


Running the simulation
-----------------------

From the directory containing ``mhd_shock_1d.py``, run:

.. code-block:: bash

    python3 -Ou mhd_shock_1d.py

The ``-O`` flag disables Python assertion overhead and ``-u`` forces unbuffered
output so progress messages appear in real time. The run completes in a few
minutes on a modern laptop. Diagnostic output is written to
``phare_outputs/shock/`` as HDF5 files.

To use multiple MPI ranks:

.. code-block:: bash

    mpirun -n 4 python3 -Ou mhd_shock_1d.py


Analyzing results
-----------------

After the run, use ``pharesee`` to load the HDF5 output and plot the solution
profiles at the final time.

.. code-block:: python

    import matplotlib.pyplot as plt
    from pyphare.pharesee.run import Run

    run = Run("phare_outputs/shock")
    t   = 80.0   # final time

    # --- density ---
    rho = run.GetMHDrho(t, merged=True)
    x_rho = rho["rho"][1][0]
    rho_vals = rho["rho"][0](x_rho)

    # --- velocity components ---
    V  = run.GetMHDV(t, merged=True)
    x_v   = V["Vx"][1][0]
    vx_vals = V["Vx"][0](x_v)

    # --- magnetic field ---
    B  = run.GetB(t, merged=True)
    x_b   = B["By"][1][0]
    by_vals = B["By"][0](x_b)

    # --- pressure ---
    P  = run.GetMHDP(t, merged=True)
    x_p   = P["P"][1][0]
    p_vals = P["P"][0](x_p)

    # Plot
    fig, axes = plt.subplots(2, 2, figsize=(11, 7), sharex=True)

    axes[0, 0].plot(x_rho, rho_vals)
    axes[0, 0].set_ylabel(r"$\rho$")

    axes[0, 1].plot(x_v, vx_vals)
    axes[0, 1].set_ylabel(r"$v_x$")

    axes[1, 0].plot(x_b, by_vals)
    axes[1, 0].set_ylabel(r"$B_y$")

    axes[1, 1].plot(x_p, p_vals)
    axes[1, 1].set_ylabel(r"$P$")

    for ax in axes.flat:
        ax.set_xlabel("x")
        ax.grid(True, alpha=0.3)

    fig.suptitle(f"Brio–Wu shock tube,  t = {t}")
    fig.tight_layout()
    plt.show()

Each ``GetMHD*`` call returns a dictionary whose keys are the component names
(``"rho"``, ``"Vx"``, ``"Vy"``, ``"Vz"``, ``"P"``). Each value is a tuple
``(interpolant, [coordinate_array])``. Evaluating ``interpolant(x)`` returns
the field values on any coordinate array you supply.

For a full description of the ``pharesee`` data API, see
:doc:`../pharesee/get_data`.


Expected output
---------------

At :math:`t = 80` the initial discontinuity has evolved into a structured
wave pattern. Reading from left to right across the domain you should see:

1. A **fast rarefaction** fan: density and pressure drop smoothly from the
   left-state values toward intermediate values.
2. A **slow shock**: a sharper drop accompanied by a rotation in :math:`B_y`.
3. A **contact discontinuity**: a density jump with no corresponding pressure jump.
4. A **slow rarefaction** or **slow shock** (depending on the beta regime): another
   partial change in the magnetic field orientation.
5. A **fast shock**: density and pressure jump to the low right-state values.

Because this run uses first-order reconstruction and the Rusanov solver the wave
transitions are smeared over several cells. Switching to ``reconstruction="linear"``
with ``riemann="hlld"`` produces sharper profiles that compare more closely to
the semi-analytical Brio–Wu reference solution.

.. note::

   The Brio–Wu problem is designed with open (outflow) boundary conditions.
   PHARE's default boundary treatment for MHD is reflective, so for very long
   run times you may see spurious reflections from the domain edges. For this
   run at :math:`t = 80` the fastest waves have not yet reached the boundaries
   and the result is unaffected.


Next steps
----------

- Switch to ``reconstruction="linear"`` with ``limiter="minmod"`` and
  ``riemann="hlld"`` to obtain second-order accuracy and cleaner wave structures.
- Enable AMR by setting ``max_nbr_levels=2`` and compare the refined and
  unrefined solutions near the shock front.
- Extend to 2D by changing ``cells`` and ``dl`` to two-element tuples and
  adding a transverse velocity perturbation to observe oblique shock formation.

See also:

- :doc:`../usage/simulation_inputs` — full reference for all configuration blocks
- :doc:`../pharesee/get_data` — reading and manipulating simulation output
- :doc:`../theory/mhd` — theoretical background of the MHD model
