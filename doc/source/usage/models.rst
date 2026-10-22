
======
Models
======

PHARE runs one of two, independent simulation types, chosen with
``Simulation(model_options=...)``:

- **Hybrid** (default, ``model_options=["HybridModel"]``): kinetic ions,
  fluid electrons. Set up with :class:`~pyphare.pharein.MaxwellianFluidModel`.
- **MHD** (``model_options=["MHDModel"]``): a single MHD fluid. Set up with
  :class:`~pyphare.pharein.MHDModel`.

There is currently no way to combine both in the same run (e.g. MHD on
coarse levels and Hybrid on fine levels) - pick one.

Hybrid: MaxwellianFluidModel
-----------------------------

The typical way to initialize a Hybrid simulation is to parametrize each ion
population by its fluid moments (density, bulk velocity, thermal velocity),
assuming an underlying Maxwellian velocity distribution.
``MaxwellianFluidModel`` also sets the magnetic field profile.

.. autoclass:: pyphare.pharein.MaxwellianFluidModel

Below, a single population of protons:

.. code-block:: python

    import pyphare.pharein as ph

    ph.MaxwellianFluidModel(
        bx=bx,
        by=by,
        bz=bz,
        protons={"charge": 1,
                 "density": density,
                 "vbulkx": vx,
                 "vbulky": vy,
                 "vbulkz": vz,
                 "vthx": vthx,
                 "vthy": vthy,
                 "vthz": vthz,
                 "nbr_part_per_cell": 100}
    )

``density``, ``vbulk*``, ``vthx/y/z`` and ``bx/by/bz`` are each either a plain
number (for a spatially uniform value) or a function of the spatial
coordinates. For example, the density below defines two current sheets at
:math:`y=0.3L_y` and :math:`y=0.7L_y`, each of half-width 0.5:

.. code-block:: python

    import numpy as np

    def density(x, y):
        # assume `sim` is the Simulation created earlier
        L = sim.simulation_domain()[1]
        return (
            0.2
            + 1.0 / np.cosh((y - L * 0.3) / 0.5) ** 2
            + 1.0 / np.cosh((y - L * 0.7) / 0.5) ** 2
        )

and a corresponding sheared magnetic field profile:

.. code-block:: python

    def S(y, y0, l):
        return 0.5 * (1.0 + np.tanh((y - y0) / l))

    def bx(x, y):
        Ly = sim.simulation_domain()[1]
        v1, v2 = -1.0, 1.0
        return v1 + (v2 - v1) * (S(y, Ly * 0.3, 0.5) - S(y, Ly * 0.7, 0.5))

Since these are plain Python functions, they can be as complex as needed
(here ``bx`` calls the helper ``S``).

Adding a second population (e.g. a lower-density beam) is done the same way,
under a different keyword:

.. code-block:: python

    ph.MaxwellianFluidModel(
        bx=bx, by=by, bz=bz,
        protons={"charge": 1, "density": density,
                 "vbulkx": vx, "vbulky": vy, "vbulkz": vz,
                 "vthx": vthx, "vthy": vthy, "vthz": vthz,
                 "nbr_part_per_cell": 100},
        beam={"charge": 1, "density": beam_density,
              "vbulkx": vx_beam, "vbulky": vy_beam, "vbulkz": vz_beam,
              "vthx": vthx_beam, "vthy": vthy_beam, "vthz": vthz_beam,
              "nbr_part_per_cell": 500},
    )

Examples:

- ``tests/functional/harris/harris_2d.py`` - one population.
- ``tests/functional/ionIonBeam/ion_ion_beam1d.py`` - two populations.

Minimum density
^^^^^^^^^^^^^^^

The electron velocity :math:`\mathbf{V}_e = (\sum_s q_s \mathbf{F}_s - \mathbf{J})/n`
and Ohm's law divide by the ion charge density, so ``Simulation`` refuses an
initial condition whose total ion charge density, in absolute value, goes below
``min_density`` (default 0.01) on any grid node of the domain, whatever the
electron closure. This is a change of behaviour: inputs with deliberately
rarefied regions must now set it explicitly, e.g.
``Simulation(..., min_density=1e-3)``.

A Hybrid simulation also requires an ``ElectronModel`` block, see below.

Hybrid: ElectronModel
-----------------------

In a Hybrid simulation, ions are kinetic but electrons are a fluid: their
momentum equation (via a chosen pressure closure) is used to compute the
electric field, assuming quasineutrality. ``ElectronModel`` sets that
closure, and is required in every Hybrid simulation.

.. autoclass:: pyphare.pharein.ElectronModel

.. code-block:: python

    import pyphare.pharein as ph

    ph.ElectronModel(closure="isothermal", Te=0.2)

Two closures are available:

- ``closure="isothermal"`` with ``Te``: the electron temperature is constant in
  space and time, and :math:`P_e = n T_e`.
- ``closure="polytropic"`` with ``Pe``, a function of the coordinates giving the
  initial electron pressure, and ``Gamma`` (default 5/3). The electrons are
  adiabatic: :math:`s = P_e / n^\gamma` is carried by the electron flow,
  :math:`\partial_t s + \mathbf{V}_e \cdot \nabla s = 0`, and
  :math:`P_e = s\, n^\gamma`. With an initial pressure
  :math:`P_e = P_0 (n/n_0)^\gamma`, :math:`s` is uniform and stays so: this is
  the usual polytropic law. A non-uniform initial :math:`s` is transported.

.. code-block:: python

    def Pe(x):
        return 0.1 * density(x) ** (5.0 / 3.0)

    ph.ElectronModel(closure="polytropic", Pe=Pe, Gamma=5.0 / 3.0)

The transport uses a first-order upwind scheme in space and is first order in
time (each predictor/corrector stage advances :math:`s` from the start of the
step), stable only if :math:`\sum_d |V_{e,d}|\, \Delta t / \Delta x_d \le 1` at
every node. The simulation stops with an error otherwise; reduce ``time_step``.
The electron pressure can be written with
``FluidDiagnostics(quantity="electron_pressure")``. Vacuum is not supported (see
the minimum density above): where the ion density drops to zero during a run,
:math:`s` and :math:`P_e` are left unchanged at that node.

Not needed for an MHD simulation.

MHD: MHDModel
--------------

.. autoclass:: pyphare.pharein.MHDModel

.. code-block:: python

    sim = ph.Simulation(
        # ...
        model_options=["MHDModel"],
        max_mhd_level=sim_max_nbr_levels,  # every level is MHD
        gamma=5.0 / 3.0,
        reconstruction="Linear",
        limiter="VanLeer",
        riemann="Rusanov",
        mhd_timestepper="TVDRK2",
    )

    ph.MHDModel(
        density=density,
        vx=vx, vy=vy, vz=vz,
        bx=bx, by=by, bz=bz,
        p=pressure,
    )

Examples: ``tests/functional/mhd_harris/harris.py``,
``tests/functional/mhd_orszagtang/orszag_tang.py``,
``tests/functional/mhd_rotor/rotor.py``.
