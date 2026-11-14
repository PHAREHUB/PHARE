=====================
The Simulation block
=====================

**Goal:** know every option of ``ph.Simulation(...)``, what it does and what
happens if you get it wrong.

``Simulation`` is the first thing a PHARE script declares. It sets how long the
run lasts, the size and resolution of the box, adaptive mesh refinement (AMR),
and where the output goes. Every other block depends on it.

.. contents:: On this page
   :local:
   :depth: 1


The blocks of a script
----------------------

A PHARE simulation is configured entirely from a Python script, by declaring a
small number of objects. Each is declared once (diagnostics may be declared
several times, once per output), and ``Simulation`` must always be declared
first.

.. code-block:: python

    # ------------ MANDATORY BLOCKS

    ph.Simulation(
        # time evolution, domain, AMR, restart and diagnostics-output parameters
    )

    ph.MaxwellianFluidModel(
        # magnetic field profile and ion population(s), for a Hybrid simulation
    )
    # -- or, for an MHD simulation --
    ph.MHDModel(
        # density, velocity, magnetic field and pressure profiles
    )

    ph.ElectronModel(
        # electron fluid closure - required for a Hybrid simulation
    )

    # ------------ END OF MANDATORY BLOCKS


    # ------------ OPTIONAL BLOCKS

    ph.LoadBalancer(
        # tune how work is redistributed across MPI ranks
    )

    ph.ElectromagDiagnostics(...)   # E and B field outputs
    ph.FluidDiagnostics(...)        # ion moment outputs
    ph.ParticleDiagnostics(...)     # particle outputs
    ph.MHDDiagnostics(...)          # MHD fluid outputs
    ph.InfoDiagnostics(...)         # run-monitoring outputs

    # ------------ END OF OPTIONAL BLOCKS

See :doc:`hybrid_inputs` for ``MaxwellianFluidModel``/``ElectronModel``,
:doc:`mhd_inputs` for ``MHDModel``, :doc:`diagnostics` for the diagnostics
blocks, and :doc:`load_balancing` for ``LoadBalancer``. This page covers the
``Simulation`` block itself.


Minimal example
---------------

.. code-block:: python

    import pyphare.pharein as ph

    ph.Simulation(
        cells=100,           # 100 cells ...
        dl=0.5,              # ... of size 0.5: the box is 50 long, in 1D
        time_step=0.01,
        time_step_nbr=1000,  # 1000 steps: the run stops at t = 10
    )

Everything else has a default. The number of space dimensions is the length of
``cells`` / ``dl`` / ``domain_size``: a number or a one-item list is 1D, two
numbers 2D, three numbers 3D.


Domain and time, worked example
-------------------------------

The domain can be given as any two of ``cells`` (number of mesh cells), ``dl``
(mesh spacing) and ``domain_size`` (physical size); the third is derived.
Likewise for time, any two of ``time_step``, ``time_step_nbr`` and
``final_time``:

.. code-block:: python

    import pyphare.pharein as ph

    ph.Simulation(
        cells=(100, 100),
        dl=(0.2, 0.2),          # domain_size is derived: (20, 20)
        time_step=0.001,
        time_step_nbr=1000,     # final_time is derived: ~1.0
        refinement="tagging",
        max_nbr_levels=3,
    )

.. plot::
   :include-source: False
   :caption: A 2D domain given by ``cells=(12, 8)`` and ``dl=(0.5, 0.5)``: its
      size ``domain_size`` is derived, here (6, 4). Dashed boxes are the
      :term:`patches <patch>` of a finer :term:`level` placed on top of it (see
      :doc:`/user_guide/amr`).

   import matplotlib.pyplot as plt
   from matplotlib.patches import Rectangle

   plt.rcParams["font.size"] = 11
   fig, ax = plt.subplots(figsize=(6.5, 4.2))
   nx, ny = 12, 8
   for i in range(nx + 1):
       ax.plot([i, i], [0, ny], color="0.8", lw=0.8)
   for j in range(ny + 1):
       ax.plot([0, nx], [j, j], color="0.8", lw=0.8)
   ax.add_patch(Rectangle((0, 0), nx, ny, fill=False, ec="0.2", lw=2))
   for (x0, y0, w, h) in [(2, 2, 4, 4), (7, 3, 3, 3)]:
       ax.add_patch(Rectangle((x0, y0), w, h, fill=False, ec="#c0392b",
                              lw=1.8, ls="--"))
   ax.text(4, 4, "patch", color="#c0392b", ha="center", va="center")
   ax.text(8.5, 4.5, "patch", color="#c0392b", ha="center", va="center")
   # one cell, with dl
   ax.annotate("", xy=(1, -0.6), xytext=(0, -0.6),
               arrowprops=dict(arrowstyle="<->", color="#2a7ab9"))
   ax.text(0.5, -1.3, "dl[0] = 0.5", color="#2a7ab9", ha="center")
   # whole domain
   ax.annotate("", xy=(nx, -2.2), xytext=(0, -2.2),
               arrowprops=dict(arrowstyle="<->", color="0.2"))
   ax.text(nx / 2, -2.9, "domain_size[0] = cells[0] * dl[0] = 12 * 0.5 = 6",
           ha="center")
   ax.annotate("", xy=(-0.7, ny), xytext=(-0.7, 0),
               arrowprops=dict(arrowstyle="<->", color="0.2"))
   ax.text(-0.9, ny / 2, "cells[1] = 8", rotation=90, ha="right", va="center")
   ax.text(nx / 2, ny + 0.4, "cells[0] = 12", ha="center")
   ax.set_xlim(-2.2, nx + 0.5)
   ax.set_ylim(-3.4, ny + 1.2)
   ax.set_aspect("equal")
   ax.axis("off")
   fig.tight_layout()

The derived values can be read back from the ``Simulation`` object, for
example to build profiles: ``sim.simulation_domain()`` gives the box size,
``sim.time_stepper.final_time`` and ``sim.time_stepper.time_step`` the time
settings.

.. note::

   When ``final_time`` and ``time_step`` are both given, ``time_step`` is
   kept exactly as given and ``final_time`` is treated as a target: the
   actual number of steps is rounded to the nearest integer, so the run may
   stop slightly before or after the requested ``final_time``. For example
   ``time_step=0.03, final_time=1.0`` gives 33 steps and stops at 0.99.

.. note::

   With ``domain_size`` and ``dl``, the number of cells is
   ``int(domain_size / dl)`` and ``dl`` is then adjusted so the box has
   exactly the size you asked: ``domain_size=10.1, dl=0.2`` gives 50 cells of
   size 0.202.


Options
-------

Each table lists the option, what it means and its default. ``cells``, ``dl``
and ``domain_size`` need one value per direction in 2D and 3D (a single number
means 1D). The patch-size options and ``nesting_buffer`` also accept a single
number, used for every direction.

Domain
~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 22 58 20

   * - Option
     - Meaning
     - Default
   * - ``cells``
     - Number of cells per direction (at least 10 in each direction).
     - none
   * - ``dl``
     - Cell size per direction, in code units.
     - none
   * - ``domain_size``
     - Box size per direction, in code units.
     - none

Give exactly two of the three.

Time
~~~~

.. list-table::
   :header-rows: 1
   :widths: 22 58 20

   * - Option
     - Meaning
     - Default
   * - ``time_step``
     - Time step :math:`dt`. A number, or a dictionary (see below).
     - none
   * - ``time_step_nbr``
     - Number of time steps.
     - none
   * - ``final_time``
     - Time at which the run stops.
     - none

Give exactly two of the three. ``time_step`` can also be a dictionary:

* ``{"mode": "constant", "value": 0.01}`` is the same as ``time_step=0.01``.
* ``{"mode": "adaptive", "cfl_wave": 0.5, "cfl_diffusive": 0.5}`` recomputes
  the time step at every step from a stability limit (:term:`CFL <CFL condition>`). It needs
  ``final_time`` and cannot be combined with ``time_step_nbr``.
  ``cfl_diffusive`` defaults to ``cfl_wave``. Both must be positive.
  Diagnostics then need explicit ``write_timestamps`` (there is no fixed
  step to build a default output grid from).

.. warning::

   Adaptive time steps work with the MHD solver only. A hybrid run with
   ``{"mode": "adaptive", ...}`` passes the Python checks, then stops at the
   start of the run with ``Adaptive time-stepping not ready in hybrid.``
   (written to the log file, the terminal only shows an MPI abort). Use a
   fixed ``time_step`` for hybrid runs.

Hybrid numerics
~~~~~~~~~~~~~~~

These apply to hybrid runs (the default).

.. list-table::
   :header-rows: 1
   :widths: 22 58 20

   * - Option
     - Meaning
     - Default
   * - ``interp_order``
     - Order of the particle shape function (:term:`interpolation order`): 1, 2
       or 3. Higher order is smoother and costlier.
     - 1
   * - ``refined_particle_nbr``
     - How many child particles a particle is split into when it enters a
       finer level (:term:`refined particle number`). Only a few values are
       valid, depending on the dimension and ``interp_order``, see below.
     - first valid value (1D: 2, 2D: 4, 3D: 6)
   * - ``particle_pusher``
     - Particle pusher (:term:`pusher`). Only ``"modified_boris"`` exists.
     - ``"modified_boris"``
   * - ``layout``
     - Where each quantity lives in the cell (:term:`Yee layout`). Only
       ``"yee"`` exists.
     - ``"yee"``
   * - ``resistivity``, ``hyper_resistivity``, ``hyper_mode``
     - Terms of Ohm's law, see :doc:`electrons`. ``hyper_mode`` is
       ``"constant"`` or ``"spatial"`` (``hyper_mode`` also applies to the
       MHD ``nu``, see :doc:`mhd_inputs`).
     - ``0.0``, ``0.0001``, ``"constant"``

Valid ``refined_particle_nbr`` values (the first of each list is the default):

.. list-table::
   :header-rows: 1
   :widths: 15 28 28 29

   * - ``interp_order``
     - 1D
     - 2D
     - 3D
   * - 1
     - 2, 3
     - 4, 5, 8, 9
     - 6, 7, 8, 9, 12 to 15, 18 to 21, 26, 27
   * - 2
     - 2, 3, 4
     - 4, 5, 8, 9, 16
     - same as order 1, and 64
   * - 3
     - 2, 3, 4, 5
     - 4, 5, 8, 9, 25
     - same as order 1, and 125

.. important::

   Python accepts all these values, but a run only works if the matching
   compiled module exists (see :doc:`/getting_started/build`). A build with
   the default combinations file contains every 1D and 2D value above, but in
   3D only 6, 12 and 27. Asking for another one stops the run at the start
   with ``This run needs the build permutation 3,1,7``, followed by the list
   of what your build has.

Patches and AMR
~~~~~~~~~~~~~~~

AMR is explained on :doc:`/user_guide/amr`. Here are only the options.

.. list-table::
   :header-rows: 1
   :widths: 22 58 20

   * - Option
     - Meaning
     - Default
   * - ``refinement``
     - How the refined regions are chosen: ``"boxes"`` (you give them, they
       never move) or ``"tagging"`` (PHARE finds them from the solution at run
       time, see :term:`tagging`). Exact lower-case spelling: other values are
       not checked.
     - ``"boxes"``
   * - ``refinement_boxes``
     - With ``"boxes"``: the fixed refined regions, per :term:`level`. The
       number of levels is deduced from it.
     - none (no refinement)
   * - ``max_nbr_levels``
     - With ``"tagging"`` only: the number of levels, level 0 included.
       Passing it with ``"boxes"`` is an error.
     - 1, with a warning (no refinement)
   * - ``tagging_threshold``
     - With ``"tagging"``: a cell is tagged when the criterion exceeds this.
       Lower means more refined cells.
     - 0.1
   * - ``tag_buffer``
     - With ``"tagging"``: number of cells added around each tagged cell
       before the :term:`patches <patch>` are built.
     - 1
   * - ``clustering``
     - How tagged cells are grouped in patches: ``"tile"`` or ``"berger"``.
     - ``"tile"``
   * - ``smallest_patch_size``
     - Smallest patch, in cells, per direction. Must be at least 3 and at most
       ``cells``. With ``interp_order=1`` the values 4 and 5 are refused (a
       known :term:`SAMRAI` bug, SAMRAI issue 311).
     - 3
   * - ``largest_patch_size``
     - Largest patch, in cells, per direction. Not smaller than the smallest
       and not larger than ``cells``.
     - no limit
   * - ``nesting_buffer``
     - Gap, in coarse cells, between the border of a level and a patch of the
       next finer level. At least 1 and at most half of ``smallest_patch_size``.
       Ignored (with a warning) if there is no refinement.
     - 1
   * - ``refinement_order``
     - Order of the operator that fills fine-level fields from the coarse
       level. Only 2 (linear) exists.
     - 2

Models
~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 22 58 20

   * - Option
     - Meaning
     - Default
   * - ``model_options``
     - ``"HybridModel"`` or ``"MHDModel"``, see :doc:`/user_guide/choosing_a_model`.
     - ``"HybridModel"``
   * - ``max_mhd_level``
     - Meant for a future run mixing MHD and hybrid levels, which does not
       work yet: leave it at its default.
     - all levels in an MHD run, 0 in a hybrid run
   * - ``gamma``, ``eta``, ``nu``, ``hall``, ``reconstruction``, ``limiter``,
       ``riemann``, ``mhd_timestepper``
     - MHD solver settings, all described on :doc:`/user_guide/mhd_inputs`.
       ``mhd_timestepper`` is mandatory for MHD.
     - see that page

.. warning::

   Running a hybrid model and an MHD model in the same run is not supported
   yet. Python accepts ``model_options=["HybridModel", "MHDModel"]``, but do
   not use it.

Output and restarts
~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 22 58 20

   * - Option
     - Meaning
     - Default
   * - ``diag_options``
     - Where and in which format diagnostics are written, see
       :doc:`/user_guide/diagnostics`.
     - none
   * - ``restart_options``
     - Checkpoints and restarting, see :doc:`/user_guide/restarts`.
       ``mode`` is mandatory when given.
     - none (no restart files)

Miscellaneous
~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 22 58 20

   * - Option
     - Meaning
     - Default
   * - ``domain_boundaries``
     - Boundary condition per box side, see
       :doc:`/user_guide/boundary_conditions`. Only MHD supports
       non-periodic sides; the box is periodic otherwise.
     - periodic
   * - ``strict``
     - Turns some warnings into errors. For instance, a model profile that is
       not periodic in a periodic box only prints
       ``Warning: Simulation is periodic but some functions are not``; with
       ``strict=True`` it stops with ``Simulation is not periodic``.
     - ``False``
   * - ``dry_run``
     - Builds and checks everything but does not start or advance the run.
       The environment variable ``PHARE_DRY_RUN=1`` has the same effect.
     - ``False``
   * - ``write_reports``
     - Writes run-environment reports in ``.phare/node_info/``: Python version
       and git hash (rank 0), and host name, memory and CPU information (one
       rank per machine).
     - ``False``
   * - ``path``
     - Stored in the simulation object; it does not choose where files are
       written (use ``diag_options`` and ``restart_options``).
     - ``"./"``
   * - ``description``
     - Free text stored in the simulation object. It is not written to the
       output files at this version.
     - none


Common mistakes
---------------

* A typo in an option name: ``interp_ordr=2`` stops with
  ``Error: invalid arguments - interp_ordr``. Fix the name.
* All three of ``cells``, ``dl`` and ``domain_size``: stops with
  ``Error: Specify either 'domain_size' and 'dl' or 'cells' and 'dl' or
  'domain_size' and 'dl'`` (the last pair repeats the first; it should read
  ``'cells' and 'domain_size'``). Remove one of the three. The same message
  appears if you give only one of them. Giving all three time options stops
  with a similar (correct) message.
* A single number mixed with tuples, e.g. ``cells=(50, 40), dl=0.5`` in 2D:
  ``Error: 'cells' and 'dl' must have the same dimension ((50, 40) and 0.5)``.
  Write ``dl=(0.5, 0.5)``.
* Fewer than 10 cells: ``Error: number of cells in non-invariant directions
  should be >= 10``.
* ``smallest_patch_size`` bigger than the box: ``Error - smallest_patch_size
  should be less than nbr of cells in all directions``.
* ``smallest_patch_size=4`` with ``interp_order=1``: ``Error:
  smallest_patch_size cannot be 4 or 5 in interp 1 because of
  https://github.com/llnl/SAMRAI/issues/311``. Use 3 or 6 and more, or
  ``interp_order=2``.
* ``refined_particle_nbr=7`` in 1D: ``Invalid split particle number, valid
  values for dim(1) interp(1) include [2, 3]``.
* A valid but not compiled ``refined_particle_nbr`` (3D, value 7): the script
  starts, then ``This run needs the build permutation 3,1,7``. Rebuild with it
  added to your combinations file, or choose a built value.
* ``max_nbr_levels`` together with the default ``refinement="boxes"``:
  ``Error: invalid arguments - max_nbr_levels``. Use ``refinement="tagging"``,
  or give ``refinement_boxes``.
* ``refinement="tagging"`` without ``max_nbr_levels``: only the warning
  ``WARNING, 'max_nbr_levels' is not set, defaulting to 1``, and no
  refinement at all.
* A typo in ``refinement`` (e.g. ``"Tagging"``) is not caught: with
  ``max_nbr_levels`` above 1 the run crashes at the start with a segmentation
  fault; without it, there is no refinement. Write exactly ``"tagging"`` or
  ``"boxes"``.
* Non-periodic ``domain_boundaries`` with a hybrid model:
  ``non-periodic domain boundaries are only supported by the MHDModel``.

.. seealso:: :doc:`/user_guide/amr`, :doc:`/user_guide/restarts`,
   :doc:`/user_guide/diagnostics`, :doc:`/getting_started/first_simulation`.


Reference
---------

.. autoclass:: pyphare.pharein.Simulation
