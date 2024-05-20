

=================
Simulation inputs
=================


Python initialization blocks
----------------------------

PHARE is initialized from Python objects.
These objects are declared by users in a script to define simulation
parameters, and to describe the initial conditionn.

The following shows the various objects to declare in a pseudo-script:

.. code-block:: python


    # ------------ MANDATORY BLOCKS

    Simulation(
        # some parameters
        # configuring numerical and AMR
        # parameters
    )


    MawellianFluidModel(
       # some parameters
       # configuring the magnetic field profile
       # and ion initial condition
       # as fluid moment profiles
       # ion particles are assumed to follow
       # locally Maxwellian distributions with these moments
       )

    # ------------ END OF MANDATORY BLOCKS


    # ------------ OPTIONAL BLOCKS

    ElectronModel(
        # configures electron fluid properties
    )



    ElectromagDiagnostics(
       # parameters configuring outputs
       # of E and B
    )


    FluidDiagnostics(
        # parameters configuring ion moment outputs
    )


    ParticleDiagnostics(
        # some parameters configuring particle outputs
    )


    # ------------ END OF OPTIONAL BLOCKS





The Simulation block
--------------------


The Simulation block is used to set general parameters of the simulation
like the integration time, domain size, interpolation order, adaptive meshing, 
restart and diagnostcs options.
The ``Simulation`` must be the first block defined in an input script

.. autoclass:: pyphare.pharein.Simulation


Magnetic field and ions
-----------------------

The typical way to initialize a Hybrid-PIC model is to parametrize the
different ion populations with fluid moments, assuming an underlying 
velocity distributotion function.
The `MaxwellianFluidModel` block allows just that, assuming a Maxwellian 
distribution for each population. 
It also allows to set the magnetic field profile.


Below is an example of a `MaxwellianFluidModel` block for which a single
population of protons is initialized.

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


In the example above, `bx`, `by`, `bz` are the components of the magnetic field,
`density` is the density of the protons, `vx`, `vy`, `vz` are the bulk velocities
and `vthx`, `vthy`, `vthz` are the thermal velocities of the protons.
All these parameters are (previously defined) functions of the spatial coordinates.


As an example, the function below defines a density profile.
This function will be called by the C++ code to load the population to which
this density is assigned.
In this example, we create two current sheets at y = 0.3 and y = 0.7 of the
simulation domain along the y direction, of half-width 0.5.

.. code-block:: python

    # assume the simulation is created and accessed by the variable `sim`

    import numpy as np

    def density(x, y):

        # L[1] is the length of the simulation domain in the y direction
        L = sim.simulation_domain()[1]
        return (
            0.2
            + 1.0 / np.cosh((y - L * 0.3) / 0.5) ** 2
            + 1.0 / np.cosh((y - L * 0.7) / 0.5) ** 2
        )


The magnetic field function could be defined as follows:

.. code-block:: python

    def S(y, y0, l):
        return 0.5 * (1.0 + np.tanh((y - y0) / l))

    def bx(x, y):
        Ly = sim.simulation_domain()[1]
        v1 = -1.0
        v2 =  1.0
        return (
            v1
            + (v2 - v1) * (S(y, Ly * 0.3, 0.5) - S(y, Ly * 0.7, 0.5))
        )

Note here that the function `bx` uses the other function `S`.
This underlines the great power that comes with initializing the simulation
with Python since the initialization script can be as complex as needed.



Adding a population is simple. In the example below a beam proton poulation
is added with more particles per cells and different moments.:


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
                 "nbr_part_per_cell": 100},

        beam={"charge": 1,
                 "density": beam_density,
                 "vbulkx": vx_beam,
                 "vbulky": vy_beam,
                 "vbulkz": vz_beam,
                 "vthx": vthx_beam,
                 "vthy": vthy_beam,
                 "vthz": vthz_beam,
                 "nbr_part_per_cell": 500}
    )


Examples
^^^^^^^^
- See tests/functional/harris/harris_2d.py for a complete example with one population.
- See tests/functional/ionIon/beam_ions.py for a complete example with two populations.

Details
^^^^^^^
- See the :doc:`maxwellian_fluid_model` class for more details.



Electron model
--------------
The `ElectronModel` block is used to set the electron fluid properties.
Below is an example of an `ElectronModel` block:

.. code-block:: python

    from pyphare.pharein import ElectronModel

    ElectronModel(closure="isothermal", Te=0.2)


For now, the only closure available is the isothermal closure, the given
temperatureis thus constantin space and time.


Diagnostics
-----------

Diagnostics blocks are used to set parameters of the different diagnostics
output PHARE produces. There are three types of diagnostics blocks:

- `ElectromagDiagnostics` for electromagnetic field outputs
- `FluidDiagnostics` for ion moment outputs
- `ParticleDiagnostics` for particle outputs



Electromagnetic Diagnostics
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Electromagnetic diagnostics are used to write either the electric or magnetic
field to disk. Below is an example of an `ElectromagDiagnostics` block.

.. code-block:: python

    from pyphare.pharein import ElectromagDiagnostics

    time_step_nbr = 1000
    time_step = 0.001
    final_time = time_step * time_step_nbr
    dt = 10 * time_step
    nt = final_time / dt + 1
    timestamps = dt * np.arange(nt)

    ElectromagDiagnostics(quantity="E", write_timestamps=timestamps)
    ElectromagDiagnostics(quantity="B", write_timestamps=timestamps)

In this example, the electric and magnetic field will be written to disk every
10 time steps. Note that although here the timestamps are equally spaced, they
do not have to be. The `write_timestamps` parameter can be any list of times.
Times are provided in simulation units.

Fluid Diagnostics
^^^^^^^^^^^^^^^^^

Fluid quantities represent the moments of the ion populations.
Quantities available are of two kinds:

Total ion quantities:
"""""""""""""""""""""

These diagnostics are properties of the ion populations as a whole.

- density: represents the total particle density of the ions
- mass density: represents the total mass density of the ions
- bulk velocity: represents the total bulk velocity of the ions
- pressure_tensory: represent the total pressure tensor of the ions

Ion population quantities:
""""""""""""""""""""""""""

These diagnostics are properties of each ion population.
The name of the population must be provided.

- flux: represents the flux of a given population
- momentum_tensor: represents the momentum tensor of a given population

Example
"""""""

.. code-block:: python

    from pyphare.pharein import FluidDiagnostics
    from pyphare.pharein import MaxwellianFluidModel

    # assume bx, by, ... are defined
    # ...
    MaxwellianFluidModel(
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

    time_step_nbr = 1000
    time_step = 0.001
    final_time = time_step * time_step_nbr
    dt = 10 * time_step
    nt = final_time / dt + 1
    timestamps = dt * np.arange(nt)

    FluidDiagnostics(
        quantity="density",
        write_timestamps=timestamps
    )

    FluidDiagnostics(
        quantity="flux",
        write_timestamps=timestamps,
        population_name="protons"
    )



Particle Diagnostics
^^^^^^^^^^^^^^^^^^^^

These diagnostics are used to write particle data to disk.
They are typically much heavier than any other diagnostics.

The block below declares a particle diagnostics so that "protons" are
written to disk at the given timestamps.
Note the `quantity="domain"` parameter. This is used to write all the particles
living within the interior of the simulation patches.

.. code-block:: python

    from pyphare.pharein import ParticleDiagnostics

    time_step_nbr = 1000
    time_step = 0.001
    final_time = time_step * time_step_nbr
    dt = 10 * time_step
    nt = final_time / dt + 1
    timestamps = dt * np.arange(nt)

    ph.ParticleDiagnostics(
            quantity="domain",
            population_name="protons",
            write_timestamps=timestamps
        )







