

=================
Simulation inputs
=================


Python script structure
-----------------------
PHARE takes python scripts as inputs. They consists in declaring various
blocks as follows.

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


The Simulation class  is used to set general parameters to the simulation
like the integration time, domain size, interpolation order, or adaptive meshing.
The ``Simulation`` must be the first block defined in an input script

.. autoclass:: pyphare.pharein.Simulation



Magnetic field and ions
-----------------------




Electron model
--------------



Diagnostics
-----------


Electromagnetic Diagnostics
^^^^^^^^^^^^^^^^^^^^^^^^^^^





Moment Diagnostics
^^^^^^^^^^^^^^^^^^




Particle Diagnostics
^^^^^^^^^^^^^^^^^^^^




Meta-data Diagnostics
^^^^^^^^^^^^^^^^^^^^^




