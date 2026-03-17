
=================
Simulation inputs
=================

PHARE is initialized from Python objects defined in a script. These objects
describe the simulation parameters, physics models, and diagnostics outputs.

.. code-block:: python

    import pyphare.pharein as ph

    # ------------ MANDATORY BLOCKS

    ph.Simulation(
        # numerical, domain, AMR, and physics parameters
    )

    ph.MaxwellianFluidModel(       # Hybrid PIC model
        # magnetic field and ion population initial conditions
    )
    # --- OR ---
    ph.MHDModel(                   # MHD model (requires model_options=["MHDModel"])
        # density, velocity, magnetic field, pressure
    )

    # ------------ OPTIONAL BLOCKS

    ph.ElectronModel(closure="isothermal", Te=0.12)   # Hybrid PIC only

    ph.ElectromagDiagnostics(quantity="B", write_timestamps=timestamps)
    ph.FluidDiagnostics(quantity="density", population_name="protons",
                        write_timestamps=timestamps)
    ph.ParticleDiagnostics(quantity="domain", population_name="protons",
                           write_timestamps=timestamps)
    ph.MHDDiagnostics(quantity="rho", write_timestamps=timestamps)

    ph.LoadBalancer(auto=True)


The Simulation block
--------------------

The ``Simulation`` block must be the first block in every input script.
It sets the domain, time stepping, numerics, AMR, physics, and output parameters.

Domain and mesh
^^^^^^^^^^^^^^^

The simulation dimensionality (1D, 2D, or 3D) is deduced from the length of
these parameters. You must specify exactly two of the three parameters
``cells``, ``dl``, and ``domain_size``; the third is computed automatically.

.. list-table::
   :header-rows: 1
   :widths: 20 15 15 50

   * - Name
     - Type
     - Default
     - Description
   * - ``cells``
     - int or tuple
     - (computed)
     - Number of cells per direction, e.g. ``100`` (1D) or ``(100, 200)`` (2D). Minimum 10 per non-invariant direction.
   * - ``dl``
     - float or tuple
     - (computed)
     - Grid spacing per direction.
   * - ``domain_size``
     - float or tuple
     - (computed)
     - Physical size of the domain per direction.
   * - ``boundary_types``
     - str or tuple
     - ``"periodic"``
     - Boundary condition type per direction. Currently only ``"periodic"`` is supported.

.. code-block:: python

    # These three are equivalent for a 1D domain of size 10 with 100 cells:
    ph.Simulation(cells=100, dl=0.1, ...)
    ph.Simulation(cells=100, domain_size=10., ...)
    ph.Simulation(dl=0.1, domain_size=10., ...)


Time stepping
^^^^^^^^^^^^^

Specify exactly two of the following three parameters; the third is computed.

.. list-table::
   :header-rows: 1
   :widths: 20 15 15 50

   * - Name
     - Type
     - Default
     - Description
   * - ``time_step``
     - float
     - (computed)
     - Duration of one time step.
   * - ``time_step_nbr``
     - int
     - (computed)
     - Total number of time steps.
   * - ``final_time``
     - float
     - (computed)
     - Simulation end time.

.. code-block:: python

    # All equivalent for 1000 steps of dt=0.001:
    ph.Simulation(time_step=0.001, time_step_nbr=1000, ...)
    ph.Simulation(time_step=0.001, final_time=1.0, ...)
    ph.Simulation(time_step_nbr=1000, final_time=1.0, ...)


Numerics
^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 20 15 15 50

   * - Name
     - Type
     - Default
     - Description
   * - ``interp_order``
     - int
     - ``1``
     - Particle b-spline interpolation order. Valid values: 1, 2, 3.
   * - ``particle_pusher``
     - str
     - ``"modified_boris"``
     - Particle pushing algorithm. Currently only ``"modified_boris"`` is available.
   * - ``layout``
     - str
     - ``"yee"``
     - Grid layout for electromagnetic fields. Currently only ``"yee"`` is supported.
   * - ``refined_particle_nbr``
     - int
     - (depends on dim/interp)
     - Number of refined particles per coarse particle during splitting. Valid values depend on dimension and interpolation order (see note below).

.. note::

   The default ``refined_particle_nbr`` is the first valid value for the given
   dimension and interpolation order. For example, in 1D with ``interp_order=1``
   the default is 2 and valid values are [2, 3]. In 2D with ``interp_order=1``
   the default is 4 and valid values are [4, 5, 8, 9].


AMR
^^^

For details on adaptive mesh refinement in PHARE, see :doc:`../theory/amr`.

.. list-table::
   :header-rows: 1
   :widths: 22 15 15 48

   * - Name
     - Type
     - Default
     - Description
   * - ``max_nbr_levels``
     - int
     - 1 (or inferred)
     - Maximum number of AMR levels (1 = no refinement). Inferred from ``refinement_boxes`` when using box-based refinement.
   * - ``refinement``
     - str
     - ``"boxes"``
     - Refinement strategy: ``"boxes"`` (static, from ``refinement_boxes``) or ``"tagging"`` (dynamic, based on ``tagging_threshold``).
   * - ``refinement_boxes``
     - dict or None
     - ``None``
     - Static refinement regions. Format: ``{"L0": {"B0": [(lo,), (up,)], ...}, ...}`` for 1D, with tuples extended for 2D/3D. Can also use ``Box`` objects.
   * - ``tagging_threshold``
     - float
     - ``0.1``
     - Threshold for dynamic tagging refinement. Only used when ``refinement="tagging"``.
   * - ``smallest_patch_size``
     - int or tuple
     - (auto)
     - Minimum number of cells per patch in each direction. Auto-computed from ``interp_order`` if not set (typically 6 or 9 depending on order, due to SAMRAI constraints).
   * - ``largest_patch_size``
     - int or tuple
     - ``None``
     - Maximum number of cells per patch in each direction. ``None`` means no upper limit.
   * - ``nesting_buffer``
     - int
     - ``0``
     - Minimum gap in coarse cells between a level boundary and any refined patch boundary.
   * - ``tag_buffer``
     - int
     - ``1``
     - Number of cells by which tagged cells are buffered before clustering. Larger values widen refined regions around tagged cells.
   * - ``clustering``
     - str
     - ``"tile"``
     - Clustering algorithm for AMR: ``"tile"`` (wider patches, better scalability) or ``"berger"`` (tighter patches).

.. code-block:: python

    # Static box-based refinement
    ph.Simulation(
        cells=100, dl=0.1,
        time_step=0.001, time_step_nbr=1000,
        refinement="boxes",
        refinement_boxes={"L0": {"B0": [(25,), (74,)]}},
    )

    # Dynamic tagging-based refinement
    ph.Simulation(
        cells=100, dl=0.1,
        time_step=0.001, time_step_nbr=1000,
        refinement="tagging",
        max_nbr_levels=3,
        tagging_threshold=0.1,
    )


Physics (Hybrid PIC)
^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 22 15 15 48

   * - Name
     - Type
     - Default
     - Description
   * - ``resistivity``
     - float
     - ``0.0``
     - Resistivity value for the generalized Ohm's law.
   * - ``hyper_resistivity``
     - float
     - ``0.0001``
     - Hyper-resistivity value. Used to control numerical noise at small scales.
   * - ``hyper_mode``
     - str
     - ``"constant"``
     - How hyper-resistivity is applied: ``"constant"`` (uniform) or ``"spatial"`` (spatially varying).

.. warning::

   The default ``hyper_resistivity`` is ``0.0001``, not zero. This provides
   baseline numerical diffusion. Set it explicitly to ``0.0`` if you want no
   hyper-resistivity.


MHD-specific parameters
^^^^^^^^^^^^^^^^^^^^^^^^

These parameters are relevant when running MHD simulations (``model_options=["MHDModel"]``).
See :doc:`../theory/mhd` for theory details.

.. list-table::
   :header-rows: 1
   :widths: 22 15 15 48

   * - Name
     - Type
     - Default
     - Description
   * - ``model_options``
     - list of str
     - ``["HybridModel"]``
     - Physics model(s) to use. Set to ``["MHDModel"]`` for MHD simulations, or ``["MHDModel", "HybridModel"]`` for coupled runs.
   * - ``gamma``
     - float
     - ``5/3``
     - Adiabatic index (ratio of specific heats).
   * - ``eta``
     - float
     - ``0.0``
     - MHD resistivity.
   * - ``nu``
     - float
     - ``0.0``
     - MHD viscosity.
   * - ``hall``
     - bool
     - ``False``
     - Enable Hall term in MHD equations.
   * - ``res``
     - bool
     - ``False``
     - Enable resistive term in MHD equations.
   * - ``hyper_res``
     - bool
     - ``False``
     - Enable hyper-resistive term in MHD equations.
   * - ``reconstruction``
     - str
     - ``""``
     - Reconstruction scheme (e.g. ``"wenoz"``, ``"mp5"``).
   * - ``limiter``
     - str
     - ``""``
     - Slope limiter (e.g. for MUSCL-type reconstruction).
   * - ``riemann``
     - str
     - ``""``
     - Riemann solver (e.g. ``"rusanov"``, ``"hlld"``).
   * - ``mhd_timestepper``
     - str
     - ``""``
     - MHD time integration scheme.
   * - ``max_mhd_level``
     - int
     - ``0``
     - Maximum AMR level on which MHD is solved (must be <= ``max_nbr_levels``).

.. note::

   Recommended MHD configurations:

   - **Robustness**: ``reconstruction="wenoz"``, ``riemann="rusanov"``
   - **Less diffusion**: ``riemann="hlld"``
   - **High accuracy**: ``reconstruction="mp5"``


Output
^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 20 15 15 50

   * - Name
     - Type
     - Default
     - Description
   * - ``diag_options``
     - dict or None
     - ``None``
     - Diagnostics output configuration. Format: ``{"format": "phareh5", "options": {"dir": "outputs/", "mode": "overwrite"}}``. Supported formats: ``"phareh5"`` (default), ``"pharevtkhdf"``.
   * - ``restart_options``
     - dict
     - ``{}``
     - Checkpoint/restart configuration. Keys: ``dir`` (path), ``mode`` (``"conserve"`` or ``"overwrite"``), ``timestamps`` (list), ``elapsed_timestamps`` (list), ``restart_time`` (float or ``"auto"``), ``keep_last`` (bool).

.. code-block:: python

    ph.Simulation(
        ...,
        diag_options={
            "format": "phareh5",
            "options": {"dir": "phare_outputs", "mode": "overwrite"},
        },
        restart_options={
            "dir": "checkpoints",
            "mode": "overwrite",
            "timestamps": [0.5, 1.0],
        },
    )


Miscellaneous
^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 20 15 15 50

   * - Name
     - Type
     - Default
     - Description
   * - ``strict``
     - bool
     - ``False``
     - If ``True``, warnings (e.g. non-periodic functions) become errors.
   * - ``dry_run``
     - bool
     - ``False``
     - If ``True``, skip model validation. Also enabled by setting the environment variable ``PHARE_DRY_RUN=1``.
   * - ``description``
     - str or None
     - ``None``
     - Arbitrary description string injected into output files when feasible.
   * - ``write_reports``
     - bool
     - ``False``
     - Write per-rank performance reports.


Complete Simulation example
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    import pyphare.pharein as ph

    ph.Simulation(
        # Domain
        cells=(200, 200),
        dl=(0.2, 0.2),
        boundary_types=("periodic", "periodic"),

        # Time
        time_step=0.001,
        time_step_nbr=10000,

        # Numerics
        interp_order=1,

        # AMR
        refinement="tagging",
        max_nbr_levels=3,
        tagging_threshold=0.1,
        clustering="tile",

        # Physics
        resistivity=0.001,
        hyper_resistivity=0.001,

        # Output
        diag_options={
            "format": "phareh5",
            "options": {"dir": "phare_outputs", "mode": "overwrite"},
        },
    )


MaxwellianFluidModel
--------------------

The ``MaxwellianFluidModel`` block defines the magnetic field profile and the
initial conditions for each ion population, assuming locally Maxwellian velocity
distributions characterized by their fluid moments. This is the standard model
for Hybrid PIC simulations.

A ``Simulation`` must be declared before creating a ``MaxwellianFluidModel``.


Magnetic field
^^^^^^^^^^^^^^

The three components ``bx``, ``by``, ``bz`` are Python callables. If omitted,
they default to ``bx=1.0``, ``by=0.0``, ``bz=0.0`` (uniform field in x).

Function signatures depend on the simulation dimension:

- **1D**: ``f(x)``
- **2D**: ``f(x, y)``
- **3D**: ``f(x, y, z)``

.. warning::

   When using periodic boundary conditions (currently the only option), all
   functions must be periodic over the simulation domain. PHARE checks this
   at initialization and will print a warning (or raise an error if
   ``strict=True``) for non-periodic functions.


Ion population parameters
^^^^^^^^^^^^^^^^^^^^^^^^^

Each ion population is passed as a keyword argument whose name becomes the
population name (e.g. ``protons={...}``). Multiple populations can be declared
in a single ``MaxwellianFluidModel`` call.

.. list-table::
   :header-rows: 1
   :widths: 22 15 15 48

   * - Name
     - Type
     - Default
     - Description
   * - ``charge``
     - float
     - ``1.0``
     - Particle charge (in units of elementary charge).
   * - ``mass``
     - float
     - ``1.0``
     - Particle mass (in units of proton mass).
   * - ``density``
     - callable
     - ``1.0`` (uniform)
     - Number density profile. Same signature as the magnetic field functions.
   * - ``vbulkx``
     - callable
     - ``0.0`` (uniform)
     - Bulk velocity in x direction.
   * - ``vbulky``
     - callable
     - ``0.0`` (uniform)
     - Bulk velocity in y direction.
   * - ``vbulkz``
     - callable
     - ``0.0`` (uniform)
     - Bulk velocity in z direction.
   * - ``vthx``
     - callable
     - ``1.0`` (uniform)
     - Thermal velocity in x direction.
   * - ``vthy``
     - callable
     - ``1.0`` (uniform)
     - Thermal velocity in y direction.
   * - ``vthz``
     - callable
     - ``1.0`` (uniform)
     - Thermal velocity in z direction.
   * - ``nbr_part_per_cell``
     - int
     - ``100``
     - Number of macro-particles per cell for this population.
   * - ``init``
     - dict
     - ``{}``
     - Initialization options. Accepts ``seed`` (int or None) for reproducible particle loading.
   * - ``density_cut_off``
     - float
     - ``1e-16``
     - Minimum density below which particles are not loaded.


Examples
^^^^^^^^

**Harris current sheet (2D, single population)**

.. code-block:: python

    import numpy as np
    import pyphare.pharein as ph

    sim = ph.Simulation(
        cells=(200, 400), dl=(0.2, 0.2),
        time_step=0.001, time_step_nbr=10000,
    )

    Ly = sim.simulation_domain()[1]

    def density(x, y):
        return (
            0.2
            + 1.0 / np.cosh((y - Ly * 0.3) / 0.5) ** 2
            + 1.0 / np.cosh((y - Ly * 0.7) / 0.5) ** 2
        )

    def S(y, y0, l):
        return 0.5 * (1.0 + np.tanh((y - y0) / l))

    def bx(x, y):
        v1, v2 = -1.0, 1.0
        return v1 + (v2 - v1) * (S(y, Ly * 0.3, 0.5) - S(y, Ly * 0.7, 0.5))

    def by(x, y):
        return 0.0

    def bz(x, y):
        return 0.0

    vth = 0.15
    ph.MaxwellianFluidModel(
        bx=bx, by=by, bz=bz,
        protons={
            "charge": 1,
            "density": density,
            "vbulkx": lambda x, y: 0.0,
            "vbulky": lambda x, y: 0.0,
            "vbulkz": lambda x, y: 0.0,
            "vthx": lambda x, y: vth,
            "vthy": lambda x, y: vth,
            "vthz": lambda x, y: vth,
            "nbr_part_per_cell": 100,
        },
    )


**Two-population example (beam + background)**

.. code-block:: python

    ph.MaxwellianFluidModel(
        bx=bx, by=by, bz=bz,
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
        },
        beam={
            "charge": 1,
            "density": beam_density,
            "vbulkx": vx_beam,
            "vbulky": vy_beam,
            "vbulkz": vz_beam,
            "vthx": vthx_beam,
            "vthy": vthy_beam,
            "vthz": vthz_beam,
            "nbr_part_per_cell": 500,
        },
    )


MHDModel
--------

The ``MHDModel`` block defines initial conditions for ideal/resistive MHD
simulations. It requires ``model_options=["MHDModel"]`` in the ``Simulation``
block.

.. list-table::
   :header-rows: 1
   :widths: 15 15 15 55

   * - Name
     - Type
     - Default
     - Description
   * - ``density``
     - callable
     - ``1.0``
     - Mass density profile.
   * - ``vx``
     - callable
     - ``1.0``
     - Velocity in x direction.
   * - ``vy``
     - callable
     - ``0.0``
     - Velocity in y direction.
   * - ``vz``
     - callable
     - ``0.0``
     - Velocity in z direction.
   * - ``bx``
     - callable
     - ``1.0``
     - Magnetic field in x direction.
   * - ``by``
     - callable
     - ``0.0``
     - Magnetic field in y direction.
   * - ``bz``
     - callable
     - ``0.0``
     - Magnetic field in z direction.
   * - ``p``
     - callable
     - ``1.0``
     - Thermal pressure.

All parameters accept callables with the same signature as ``MaxwellianFluidModel``
functions (``f(x)`` in 1D, ``f(x,y)`` in 2D, ``f(x,y,z)`` in 3D). If omitted,
the listed defaults are used as uniform values.

.. note::

   Recommended MHD solver configurations:

   - **Robustness**: ``reconstruction="wenoz"``, ``riemann="rusanov"``
   - **Less numerical diffusion**: ``riemann="hlld"``
   - **High-order accuracy**: ``reconstruction="mp5"``

.. code-block:: python

    import numpy as np
    import pyphare.pharein as ph

    ph.Simulation(
        cells=400, dl=0.1,
        time_step=0.0001, time_step_nbr=10000,
        model_options=["MHDModel"],
        gamma=5.0 / 3.0,
        reconstruction="wenoz",
        riemann="rusanov",
        mhd_timestepper="euler",
    )

    ph.MHDModel(
        density=lambda x: 1.0 + 0.5 * np.sin(2 * np.pi * x / 40.0),
        vx=lambda x: 0.0,
        vy=lambda x: 0.0,
        vz=lambda x: 0.0,
        bx=lambda x: 1.0,
        by=lambda x: 0.0,
        bz=lambda x: 0.0,
        p=lambda x: 1.0,
    )


ElectronModel
-------------

The ``ElectronModel`` block sets the electron fluid closure for Hybrid PIC
simulations. It is optional but recommended.

.. list-table::
   :header-rows: 1
   :widths: 15 15 15 55

   * - Name
     - Type
     - Default
     - Description
   * - ``closure``
     - str
     - (required)
     - Closure model. Currently only ``"isothermal"`` is supported.
   * - ``Te``
     - float
     - ``0.1``
     - Electron temperature. The electron pressure is computed as Pe = n * Te.

.. code-block:: python

    ph.ElectronModel(closure="isothermal", Te=0.12)

.. note::

   The isothermal closure means the electron temperature is constant in both
   space and time. The electron pressure contribution to the electric field is
   Pe = n * Te, where n is the local ion charge density.


Diagnostics
-----------

Diagnostics blocks configure what data PHARE writes to disk during a simulation.
All diagnostics share these common parameters:

.. list-table::
   :header-rows: 1
   :widths: 22 15 15 48

   * - Name
     - Type
     - Default
     - Description
   * - ``write_timestamps``
     - array-like
     - (required)
     - Simulation times at which to write data. Does not need to be uniformly spaced. Values must be multiples of ``time_step``.
   * - ``quantity``
     - str
     - (required)
     - Which quantity to output (type-specific, see below).
   * - ``flush_every``
     - int
     - ``1``
     - Flush HDF5 buffers every N dumps. ``1`` is safest, higher values improve I/O performance. ``0`` means never flush explicitly.

.. note::

   ``write_timestamps`` is the primary mechanism for controlling when data is
   written. Timestamps are given in simulation time units and must be consistent
   with ``time_step``.


ElectromagDiagnostics
^^^^^^^^^^^^^^^^^^^^^

Writes electromagnetic field data.

- ``quantity="E"`` -- electric field (Ex, Ey, Ez)
- ``quantity="B"`` -- magnetic field (Bx, By, Bz)

.. code-block:: python

    import numpy as np

    dt = 10 * time_step
    timestamps = dt * np.arange(int(final_time / dt) + 1)

    ph.ElectromagDiagnostics(quantity="E", write_timestamps=timestamps)
    ph.ElectromagDiagnostics(quantity="B", write_timestamps=timestamps)


FluidDiagnostics
^^^^^^^^^^^^^^^^

Writes ion fluid moment data. Quantities are divided into total-ion and
per-population categories.

**Total ion quantities** (no ``population_name`` needed):

- ``"charge_density"`` -- total ion charge density
- ``"mass_density"`` -- total ion mass density
- ``"bulkVelocity"`` -- total ion bulk velocity
- ``"momentum_tensor"`` -- total ion momentum tensor

**Per-population quantities** (``population_name`` required):

- ``"density"`` -- population number density
- ``"flux"`` -- population particle flux
- ``"momentum_tensor"`` -- population momentum tensor

**Special wrapper**:

- ``"pressure_tensor"`` -- automatically creates the diagnostics needed to
  compute the pressure tensor. For total ions this creates ``mass_density``,
  ``bulkVelocity``, and ``momentum_tensor``. For a specific population it
  creates ``density``, ``flux``, and ``momentum_tensor``.

.. code-block:: python

    # Total ion bulk velocity
    ph.FluidDiagnostics(quantity="bulkVelocity", write_timestamps=timestamps)

    # Per-population density
    ph.FluidDiagnostics(
        quantity="density",
        population_name="protons",
        write_timestamps=timestamps,
    )

    # Pressure tensor shortcut (creates multiple diagnostics)
    ph.FluidDiagnostics(quantity="pressure_tensor", write_timestamps=timestamps)

.. warning::

   ``"flux"`` is only available for specific populations, not for total ions
   (use ``"bulkVelocity"`` instead). ``"density"`` always requires a
   ``population_name`` (use ``"charge_density"`` or ``"mass_density"`` for
   total ions).


ParticleDiagnostics
^^^^^^^^^^^^^^^^^^^

Writes raw particle data. These files are typically much larger than field
or fluid diagnostics. ``population_name`` is always required.

- ``"domain"`` -- all particles in patch interiors (most common)
- ``"levelGhost"`` -- particles in level ghost regions
- ``"space_box"`` -- particles within a user-specified spatial region (requires ``extent`` parameter)

.. code-block:: python

    ph.ParticleDiagnostics(
        quantity="domain",
        population_name="protons",
        write_timestamps=timestamps,
    )


MHDDiagnostics
^^^^^^^^^^^^^^

Writes MHD state variables. Only available when running MHD simulations.

- ``"rho"`` -- mass density
- ``"V"`` -- velocity
- ``"P"`` -- pressure
- ``"rhoV"`` -- momentum density
- ``"Etot"`` -- total energy density

.. code-block:: python

    ph.MHDDiagnostics(quantity="rho", write_timestamps=timestamps)
    ph.MHDDiagnostics(quantity="V", write_timestamps=timestamps)
    ph.MHDDiagnostics(quantity="P", write_timestamps=timestamps)


MetaDiagnostics and InfoDiagnostics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- ``MetaDiagnostics(quantity="tags", ...)`` -- writes AMR refinement tags.
- ``InfoDiagnostics(quantity="particle_count", ...)`` -- writes particle
  count per rank. If ``write_timestamps`` is omitted, it defaults to every
  simulation time step.


Output formats
^^^^^^^^^^^^^^

The output format is set via ``diag_options`` in the ``Simulation`` block:

- ``"phareh5"`` (default) -- native PHARE HDF5 format
- ``"pharevtkhdf"`` -- VTK HDF format for visualization in ParaView

.. code-block:: python

    ph.Simulation(
        ...,
        diag_options={"format": "pharevtkhdf", "options": {"dir": "vtk_output"}},
    )


LoadBalancer
------------

The ``LoadBalancer`` controls dynamic load balancing across MPI ranks during
the simulation. It is optional; if omitted, no rebalancing is performed.

.. list-table::
   :header-rows: 1
   :widths: 22 15 15 48

   * - Name
     - Type
     - Default
     - Description
   * - ``active``
     - bool
     - ``True``
     - Whether load balancing is enabled.
   * - ``mode``
     - str
     - ``"nppc"``
     - How load is assessed: ``"nppc"`` (particle count per rank) or ``"homogeneous"`` (cell count per rank).
   * - ``tol``
     - float
     - ``0.05``
     - Acceptable load imbalance tolerance.
   * - ``on_init``
     - bool
     - ``False``
     - Whether to rebalance at initialization.
   * - ``auto``
     - bool
     - ``False``
     - Automatic rebalancing with exponential backoff. Mutually exclusive with ``every``.
   * - ``every``
     - int or None
     - ``None``
     - Rebalance every N time steps. If ``None``, ``auto`` mode is activated.

.. note::

   If ``every`` is not set, ``auto`` mode is automatically enabled. You
   cannot set both ``every`` and ``auto=True``.

**Expert parameters** (used only with ``auto=True``):

.. list-table::
   :header-rows: 1
   :widths: 35 10 10 45

   * - Name
     - Type
     - Default
     - Description
   * - ``next_rebalance``
     - int
     - ``200``
     - Initial number of steps before first auto-rebalance check.
   * - ``next_rebalance_backoff_multiplier``
     - int
     - ``2``
     - Multiplier for exponential backoff between rebalance checks.
   * - ``max_next_rebalance``
     - int
     - ``1000``
     - Maximum interval between rebalance checks.

.. code-block:: python

    # Simple: rebalance every 500 steps
    ph.LoadBalancer(every=500, mode="nppc", tol=0.05)

    # Auto mode with defaults
    ph.LoadBalancer(auto=True)
