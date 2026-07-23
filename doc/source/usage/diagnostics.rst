
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
``write_timestamps``/``elapsed_timestamps`` (when to write) and ``quantity``
(what to write - valid values differ per block, see below).

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
need to be equally spaced, and are given in simulation time units.

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
   ``"pharevtkhdf"``) is ``diag_options["format"]``.

Particle Diagnostics
----------------------

.. autoclass:: pyphare.pharein.ParticleDiagnostics

MHD Diagnostics
-----------------

.. autoclass:: pyphare.pharein.MHDDiagnostics

Meta Diagnostics
------------------

.. autoclass:: pyphare.pharein.MetaDiagnostics

Info Diagnostics
------------------

.. autoclass:: pyphare.pharein.InfoDiagnostics
