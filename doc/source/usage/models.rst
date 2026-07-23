
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

``density``, ``vbulk*``, ``vthx/y/z`` and ``bx/by/bz`` are all functions of
the spatial coordinates. For example, the density below defines two current
sheets at :math:`y=0.3L_y` and :math:`y=0.7L_y`, each of half-width 0.5:

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
- ``tests/functional/ionIon/beam_ions.py`` - two populations.

A Hybrid simulation also requires an :doc:`ElectronModel <electrons>` block.

MHD: MHDModel
--------------

.. autoclass:: pyphare.pharein.MHDModel

.. code-block:: python

    sim = ph.Simulation(
        # ...
        model_options=["MHDModel"],
        max_mhd_level=sim_max_nbr_levels,  # every level is MHD
        gamma=5.0 / 3.0,
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
