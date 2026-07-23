
==========
Simulation
==========

Setting up a simulation
------------------------

A PHARE simulation is configured entirely from a Python script, by declaring
a small number of objects. Each is declared once (diagnostics may be declared
several times, once per output), and ``Simulation`` must always be declared
first - every other block below depends on it.

.. code-block:: python

    # ------------ MANDATORY BLOCKS

    Simulation(
        # time evolution, domain, AMR, restart and diagnostics-output parameters
    )

    MaxwellianFluidModel(
        # magnetic field profile and ion population(s), for a Hybrid simulation
    )
    # -- or, for an MHD simulation --
    MHDModel(
        # density, velocity, magnetic field and pressure profiles
    )

    ElectronModel(
        # electron fluid closure - required for a Hybrid simulation
    )

    # ------------ END OF MANDATORY BLOCKS


    # ------------ OPTIONAL BLOCKS

    LoadBalancer(
        # tune how work is redistributed across MPI ranks
    )

    ElectromagDiagnostics(...)   # E and B field outputs
    FluidDiagnostics(...)        # ion moment outputs
    ParticleDiagnostics(...)     # particle outputs
    MHDDiagnostics(...)          # MHD fluid outputs
    MetaDiagnostics(...)         # AMR bookkeeping outputs
    InfoDiagnostics(...)         # run-monitoring outputs

    # ------------ END OF OPTIONAL BLOCKS

See :doc:`models` for ``MaxwellianFluidModel``/``MHDModel``, :doc:`electrons`
for ``ElectronModel``, :doc:`diagnostics` for the diagnostics blocks, and
:doc:`load_balancing` for ``LoadBalancer``. This page covers the
``Simulation`` block itself.

The Simulation block
---------------------

``Simulation`` sets the overall numerical and physical parameters of the
run: how long it evolves and at what resolution, the physical domain and
mesh, adaptive mesh refinement, restart checkpoints, and where diagnostics
are written.

.. autoclass:: pyphare.pharein.Simulation

Domain and time, worked example
--------------------------------

The domain can be given as any two of ``cells`` (number of mesh cells),
``dl`` (mesh spacing) and ``domain_size`` (physical size) - the third is
derived. Likewise for time, any two of ``time_step``, ``time_step_nbr`` and
``final_time``:

.. code-block:: python

    from pyphare.pharein import Simulation

    Simulation(
        cells=(100, 100),
        dl=(0.2, 0.2),          # domain_size is derived: (20, 20)
        time_step=0.001,
        time_step_nbr=1000,     # final_time is derived: ~1.0
        max_nbr_levels=3,
    )

.. note::

   When ``final_time`` and ``time_step`` are both given, ``time_step`` is
   kept exactly as given and ``final_time`` is treated as a target: the
   actual number of steps is rounded to the nearest integer, so the run may
   stop slightly before or after the requested ``final_time``.

Adaptive mesh refinement
-------------------------

Refined regions can be placed in two ways, set by ``refinement``:

- ``refinement="boxes"`` (default): you give the exact boxes to refine, one
  list per level, in ``refinement_boxes``. Boxes never move.
- ``refinement="tagging"``: PHARE decides where to refine at run time, based
  on the local solution. This is what most functional examples in
  ``tests/functional/`` use.

.. code-block:: python

    # static boxes, one region on level 0
    Simulation(
        # ...
        refinement_boxes={"L0": {"B0": [(10, 10), (50, 50)]}},
    )

    # dynamic, gradient-based refinement
    Simulation(
        # ...
        refinement="tagging",
        max_nbr_levels=3,
        tagging_threshold=0.1,
        tag_buffer=2,
    )

.. warning::

   With ``interp_order=1``, ``smallest_patch_size`` cannot be 4 or 5 in any
   direction (a known upstream SAMRAI limitation). Pick a different value,
   e.g. 6 or larger, or leave ``smallest_patch_size`` at its default.

Restarts and diagnostics output
---------------------------------

See :doc:`restarts` for ``restart_options`` and :doc:`diagnostics` for
``diag_options`` and the diagnostics blocks themselves.
