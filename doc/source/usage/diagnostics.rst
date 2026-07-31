
===========
Diagnostics
===========

Diagnostics blocks tell PHARE what to write to disk, and when. Several
diagnostics blocks are available, one per family of physical quantity; each
can be declared multiple times (e.g. once per quantity, or once per
population). Declaring the same quantity (and, for per-population
diagnostics, the same population) twice does not create a duplicate: it
just extends the existing one's list of output times.

All diagnostics blocks share the parameters documented on
:class:`~pyphare.pharein.diagnostics.Diagnostics`, in particular
``quantity`` (what to write - valid values differ per block, see below) and,
controlling when to write, ``write_timestamps`` and/or
``elapsed_timestamps``. At least one of the two is required.

``write_timestamps`` is a list of simulation times, in simulation time
units, and each value must line up with the simulation's ``time_step``:

.. code-block:: python

    import numpy as np

    time_step_nbr = 1000
    time_step = 0.001
    final_time = time_step * time_step_nbr
    dt = 10 * time_step
    timestamps = dt * np.arange(final_time / dt + 1)

    ph.ElectromagDiagnostics(quantity="E", write_timestamps=timestamps)
    ph.ElectromagDiagnostics(quantity="B", write_timestamps=timestamps)

In this example E and B are written every 10 time steps. Timestamps don't
need to be equally spaced.

``elapsed_timestamps`` is the wall-clock counterpart: instead of simulation
times, it is a list of wall-clock durations (in seconds, or as
``datetime.timedelta``) elapsed since the simulation started. It is handy
when you don't know in advance which simulation time a run will reach within
a given wall-clock budget (e.g. on a cluster with a time limit), or when a
diagnostic should be dumped on a wall-clock cadence rather than a
simulation-time one:

.. code-block:: python

    import datetime

    ph.ElectromagDiagnostics(
        quantity="E",
        elapsed_timestamps=[datetime.timedelta(minutes=10 * i) for i in range(6)],
    )

This writes ``E`` roughly every 10 minutes of wall-clock time, for the first
hour of the run. Unlike ``write_timestamps``, each dump is triggered by
comparing against the actual wall-clock time as the run progresses, so the
simulation time (and time step) at which it happens depends on how fast the
run goes - it is not reproducible across runs at different speeds or on
different machines. ``write_timestamps`` and ``elapsed_timestamps`` can be
combined on the same diagnostics block; PHARE writes whenever either
condition is met.

Where diagnostics files are written, and in what format, is controlled by
``Simulation(diag_options=...)`` (see :doc:`simulation`) rather than by the
diagnostics blocks themselves.

Electromagnetic Diagnostics
----------------------------

.. autoclass:: pyphare.pharein.ElectromagDiagnostics

Fluid Diagnostics
------------------

.. autoclass:: pyphare.pharein.FluidDiagnostics

.. important::

   ``Simulation(diag_export_format=...)`` and ``diag_options["format"]`` are
   two different, unrelated settings despite the similar name:
   ``diag_export_format`` currently only accepts ``"hdf5"`` and has no other
   effect; the actual on-disk diagnostics format (``"phareh5"`` or
   ``"pharevtkhdf"``) is ``diag_options["format"]``, described below.

.. admonition:: File formats: phareh5 vs pharevtkhdf

   ``"phareh5"`` is PHARE's own `HDF5 <https://www.hdfgroup.org/solutions/hdf5/>`_-based
   format. Quantities are written with their exact centering on the mesh
   (including ghost nodes), which is what makes it the format used in most
   tests: it preserves the raw, uninterpolated data.

   ``"pharevtkhdf"`` is Kitware's
   `VTKHDF <https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html#vtkhdf-file-format>`_
   format. All quantities are interpolated to mesh nodes (point values),
   and the files are readily readable by ParaView.

Particle Diagnostics
----------------------

.. autoclass:: pyphare.pharein.ParticleDiagnostics

MHD Diagnostics
-----------------

.. autoclass:: pyphare.pharein.MHDDiagnostics

Info Diagnostics
------------------

.. autoclass:: pyphare.pharein.InfoDiagnostics
