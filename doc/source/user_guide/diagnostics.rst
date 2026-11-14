===========
Diagnostics
===========

**Goal:** tell PHARE which quantities to write to disk, and at which times, then
find the files it produced.

Minimal example
---------------

Add these lines to a simulation script (the full script is the one from
:doc:`/getting_started/first_simulation`, with a shorter run). Here the fields,
the ion densities, the ion velocity and the particles are written at three times:

.. code-block:: python

    times = [0.0, 0.5, 1.0]  # multiples of time_step

    ph.ElectromagDiagnostics(quantity="E", write_timestamps=times)
    ph.ElectromagDiagnostics(quantity="B", write_timestamps=times)
    ph.FluidDiagnostics(quantity="charge_density", write_timestamps=times)
    ph.FluidDiagnostics(quantity="bulkVelocity", write_timestamps=times)
    ph.FluidDiagnostics(quantity="density", population_name="protons", write_timestamps=times)
    ph.FluidDiagnostics(quantity="flux", population_name="protons", write_timestamps=times)
    ph.ParticleDiagnostics(quantity="domain", population_name="protons",
                           write_timestamps=[0.0, 1.0])

With ``diag_options={"format": "phareh5", "options": {"dir": "diag_run",
"mode": "overwrite"}}`` in ``Simulation(...)`` and a run to ``final_time=1.0``
(100 steps, 100 particles per cell, about 2 seconds), the folder ``diag_run``
contains:

.. code-block:: text

    EM_B.h5                    EM_E.h5
    ions_bulkVelocity.h5       ions_charge_density.h5
    ions_pop_protons_density.h5  ions_pop_protons_flux.h5
    ions_pop_protons_domain.h5

What happens
------------

Each diagnostics block names one quantity. At every time you listed, PHARE
appends that quantity (all patches, all refinement levels) to **one file per
quantity**. The file name is the quantity path with ``/`` replaced by ``_``.
You read the files back with ``pharesee`` (see :doc:`/data_analysis/reading_data`).

Declaring the same quantity (and, for per-population diagnostics, the same
population) twice does not create a duplicate: it just extends the existing
one's list of output times (only ``write_timestamps`` are merged). Each block can be declared as many times as you
like, for instance once per quantity or once per population.

Which blocks exist
------------------

.. list-table::
   :header-rows: 1
   :widths: 30 45 25

   * - Block
     - Writes
     - Written by ``pharevtkhdf``?
   * - ``ElectromagDiagnostics``
     - electric and magnetic fields
     - yes
   * - ``FluidDiagnostics``
     - ion moments (:term:`hybrid model`)
     - yes
   * - ``ParticleDiagnostics``
     - macro-particles of one population
     - no
   * - ``MHDDiagnostics``
     - :term:`MHD` fields
     - yes
   * - ``InfoDiagnostics``
     - run monitoring (particle count)
     - no
   * - ``MetaDiagnostics``
     - refinement tags
     - no

The "no" blocks are accepted, but ``pharevtkhdf`` silently writes nothing for
them. If you asked for
particles and see no file, check ``diag_options["format"]``.

All blocks take the parameters documented on
:class:`~pyphare.pharein.diagnostics.Diagnostics`: ``quantity`` (what to
write) and ``write_timestamps`` and/or ``elapsed_timestamps`` (when to write;
at least one is required), plus ``flush_every``. (``path`` is accepted but has
no effect: the folder is set by ``diag_options``, below.)

Valid quantities
----------------

Electromagnetic fields
~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 20 40 40

   * - ``quantity``
     - Meaning
     - File
   * - ``"E"``
     - electric field
     - ``EM_E.h5``
   * - ``"B"``
     - magnetic field
     - ``EM_B.h5``

Ion fluid moments (hybrid model)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For the ions as a whole, do **not** give ``population_name``:

.. list-table::
   :header-rows: 1
   :widths: 22 38 40

   * - ``quantity``
     - Meaning
     - File
   * - ``"charge_density"``
     - charge density
     - ``ions_charge_density.h5``
   * - ``"mass_density"``
     - mass density
     - ``ions_mass_density.h5``
   * - ``"bulkVelocity"``
     - bulk velocity of the ions
     - ``ions_bulkVelocity.h5``
   * - ``"charge_flux"``
     - charge flux (ion current density)
     - ``ions_charge_flux.h5``
   * - ``"momentum_tensor"``
     - momentum tensor
     - ``ions_momentum_tensor.h5``

For one population, give ``population_name="protons"`` (the name used in
:class:`~pyphare.pharein.MaxwellianFluidModel`):

.. list-table::
   :header-rows: 1
   :widths: 22 34 44

   * - ``quantity``
     - Meaning
     - File
   * - ``"density"``
     - particle density of the population
     - ``ions_pop_protons_density.h5``
   * - ``"charge_density"``
     - charge density of the population
     - ``ions_pop_protons_charge_density.h5``
   * - ``"flux"``
     - flux of the population
     - ``ions_pop_protons_flux.h5``
   * - ``"momentum_tensor"``
     - momentum tensor of the population
     - ``ions_pop_protons_momentum_tensor.h5``

``density`` and ``flux`` require ``population_name``; ``charge_flux`` is only
for the ions as a whole.

Particles
~~~~~~~~~

``ParticleDiagnostics`` always needs ``population_name``.

.. list-table::
   :header-rows: 1
   :widths: 18 46 36

   * - ``quantity``
     - Meaning
     - File
   * - ``"domain"``
     - particles inside the patches (default, use this)
     - ``ions_pop_protons_domain.h5``
   * - ``"levelGhost"``
     - particles in the ghost zone between two levels, for debugging refinement
     - ``ions_pop_protons_levelGhost.h5``

Particle files are much larger than field files: write them rarely.

MHD model
~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 20 40 40

   * - ``quantity``
     - Meaning
     - File
   * - ``"rho"``
     - mass density
     - ``mhd_rho.h5``
   * - ``"V"``
     - velocity
     - ``mhd_V.h5``
   * - ``"P"``
     - thermal pressure
     - ``mhd_P.h5``
   * - ``"rhoV"``
     - momentum density
     - ``mhd_rhoV.h5``
   * - ``"Etot"``
     - total energy density
     - ``mhd_Etot.h5``

Run monitoring
~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 20 50 30

   * - ``quantity``
     - Meaning
     - File
   * - ``"particle_count"``
     - number of particles per patch, to follow the balance between patches over the run
     - ``particle_count.h5``

``MetaDiagnostics(quantity="tags")`` writes the refinement tags (cells marked
for refinement) to ``tags.h5``.

``InfoDiagnostics`` writes at every time step if you give no
``write_timestamps`` (fixed time step only).

With ``pharevtkhdf`` the file extension is ``.vtkhdf`` instead of ``.h5``.

When to write: timestamps
-------------------------

``write_timestamps`` is a list of simulation times, in simulation time units.
Use **multiples of** ``time_step``: a diagnostic can only be written at a time
the simulation actually reaches. The rules:

* a time below the start time is dropped, with the message
  ``Warning: some timestamps below (0) are filtered``;
* a time above ``final_time`` is an error;
* the list must be in ascending order (otherwise ``RuntimeError``);
* the times need not be equally spaced.

.. code-block:: python

    import numpy as np
    import pyphare.pharein as ph

    time_step_nbr = 1000
    time_step = 0.001
    final_time = time_step * time_step_nbr
    dt = 10 * time_step
    timestamps = dt * np.arange(final_time / dt + 1)

    ph.ElectromagDiagnostics(quantity="E", write_timestamps=timestamps)
    ph.ElectromagDiagnostics(quantity="B", write_timestamps=timestamps)

In this example E and B are written every 10 time steps.

``elapsed_timestamps`` is the wall-clock counterpart: durations since the run
started (seconds or ``datetime.timedelta``). Use it on a cluster, when you do
not know which simulation time the job will reach:

.. code-block:: python

    import datetime

    ph.ElectromagDiagnostics(
        quantity="E",
        elapsed_timestamps=[datetime.timedelta(minutes=10 * i) for i in range(6)],
    )

This writes ``E`` about every 10 minutes for the first hour. The simulation
time of each dump depends on the speed of the run, so it is not reproducible.
Both lists can be given on one block: PHARE writes when either is due.

Where and how files are written: ``diag_options``
-------------------------------------------------

Where diagnostics files are written, and in what format, is controlled by
``Simulation(diag_options=...)`` (see :doc:`simulation`), not by the
diagnostics blocks themselves.

.. list-table::
   :header-rows: 1
   :widths: 28 46 26

   * - Key
     - Meaning
     - Default
   * - ``format``
     - ``"phareh5"`` or ``"pharevtkhdf"``
     - ``"phareh5"``
   * - ``options["dir"]``
     - directory for the files (created if missing, when ``format`` and
       ``dir`` are both given)
     - ``"phare_output"`` when ``diag_options`` is not given, but PHARE then
       does **not create** the folder: always give ``format`` and
       ``options["dir"]``
   * - ``options["mode"]``
     - only ``"overwrite"``: each diagnostics file PHARE writes is emptied when
       it is opened at the start of the run (other files in ``dir`` are left alone)
     - not set (see below)
   * - ``options["allow_emergency_dumps"]``
     - ``True`` or ``False``: dump everything if the run stops on an error
     - ``False``

* Without ``mode``, an existing file is opened and kept, and new data is
  added to it (useful when restarting). With ``"overwrite"`` it starts empty.
* An *emergency dump* writes all diagnostics of all refinement levels, outside
  the schedule, at the current time, when the run stops on an error in the
  particle update, in the initialisation of the hybrid levels, or in the MHD
  solver. Use it to look at the state just before a crash.
* An invalid ``format`` or ``mode``, or a non-boolean ``allow_emergency_dumps``,
  raises ``ValueError`` (checked only when ``format`` is given: always give it).

.. important::

   The file format is chosen only with ``diag_options["format"]``. The
   ``Simulation`` reference (:doc:`simulation`) also lists a ``diag_export_format`` argument,
   but ``Simulation`` does not accept it: passing it stops the script with
   ``ValueError: Error: invalid arguments - diag_export_format``.

.. admonition:: File formats: phareh5 vs pharevtkhdf

   ``"phareh5"`` is PHARE's own `HDF5 <https://www.hdfgroup.org/solutions/hdf5/>`_-based
   format. Quantities are written with their exact centering on the mesh
   (including ghost nodes), which is what makes it the format used in most
   tests: it preserves the raw, uninterpolated data. It is the format read by
   ``pharesee``'s ``Run``.

   ``"pharevtkhdf"`` is Kitware's
   `VTKHDF <https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html#vtkhdf-file-format>`_
   format. All quantities are interpolated to mesh nodes (point values),
   and the files are readily readable by ParaView. It only writes the fluid,
   electromagnetic and MHD blocks (see above).

Precision
---------

Diagnostics are written in **single precision**, in both formats, unless PHARE
was built with ``-DPHARE_DIAG_DOUBLES=1``. Single precision keeps about 7
significant digits, so very small perturbations, differences of nearly equal
fields or divergence checks can be limited by the file, not by the
simulation. See :doc:`/getting_started/build` to switch it on.

Common mistakes
---------------

* ``Error: 'Bx' not a valid electromag diagnostics : E, B`` — wrong
  ``quantity``. This is only **printed**: the script goes on and that
  diagnostic is silently not registered, so read the start of your log.
  The same holds for the next messages.
* ``Error: cannot use density without population name`` (for ``flux``:
  ``Error: missing population_name``) — ``density`` and ``flux`` need
  ``population_name="protons"``.
* ``AttributeError: 'NoneType' object has no attribute 'populations'`` — declare
  ``MaxwellianFluidModel`` before the diagnostics that use ``population_name``.
* ``RuntimeError: Error: ElectromagDiagnostics.write_timestamps is inconsistent
  with simulation.time_step`` — a time that is not a multiple of
  ``time_step``. Some such times are not caught, so check your list yourself.
* ``'charge_flux' is only available for the ions as a whole`` — remove
  ``population_name``.
* ``Error: population 'alpha' not in simulation initial model`` — the name
  must be one declared in ``MaxwellianFluidModel``.
* ``RuntimeError: Error: missing parameters - one required (non-empty):
  write_timestamps, elapsed_timestamps`` — give at least one of them.
* ``RuntimeError: Error: timestamp([5.]) cannot be greater than
  simulation.final_time(1.0))`` — a time beyond the end of the run.
* ``RuntimeError: Error: ElectromagDiagnostics.write_timestamps not in
  ascending order)`` — sort the list.
* Asked for particles with ``pharevtkhdf`` and got no file — that format
  does not write them; use ``phareh5``.
* ``mass_density`` or ``bulkVelocity`` with ``population_name`` is accepted by
  Python (and the ``FluidDiagnostics`` reference below says it works), but the
  run then stops (``Unknown Diagnostic Quantity``): use them for the ions as a
  whole.
* ``Failed to create file phare_output/ions_charge_density.h5 (... MPI function
  failed)`` — the output folder does not exist. PHARE creates it only when
  ``diag_options`` gives both ``format`` and ``options["dir"]``: give both.

Reference
---------

Electromagnetic Diagnostics
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: pyphare.pharein.ElectromagDiagnostics

Fluid Diagnostics
~~~~~~~~~~~~~~~~~

.. autoclass:: pyphare.pharein.FluidDiagnostics

Particle Diagnostics
~~~~~~~~~~~~~~~~~~~~

.. autoclass:: pyphare.pharein.ParticleDiagnostics

MHD Diagnostics
~~~~~~~~~~~~~~~

.. autoclass:: pyphare.pharein.MHDDiagnostics

Info Diagnostics
~~~~~~~~~~~~~~~~

.. autoclass:: pyphare.pharein.InfoDiagnostics

Meta Diagnostics
~~~~~~~~~~~~~~~~

.. autoclass:: pyphare.pharein.MetaDiagnostics
