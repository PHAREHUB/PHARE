=====================
Environment variables
=====================

PHARE reads a number of environment variables at runtime, on both the C++ and
Python sides.

``PHARE_LOG``
--------------

Redirects C++ ``std::cout`` output.

- ``RANK_FILES`` (default) - one log file per MPI rank, under ``.log/``
- ``DATETIME_FILES`` - one log file per rank, named with a startup datetime timestamp
- ``CLI`` - no redirection, prints to the terminal
- ``NULL`` - discarded (sent to ``/dev/null``)

Any other value is rejected. When running via ``pyphare``, the Python side defaults this to
``RANK_FILES`` if unset (and creates the ``.log/`` directory on rank 0) before the C++ simulator
starts. When running pure C++ (e.g. unit tests) without going through Python, no redirection
happens unless this variable is set.

``PHARE_SCOPE_TIMING``
-----------------------

Enables ``phlop`` scope/timing instrumentation. Truthy values: ``1`` or ``true``.

Writes one timing file per rank to ``.phare/timings/rank.<N>.txt``. Read independently on both
the C++ side (gated behind the ``phlop`` build option) and the Python side (which additionally
triggers monitoring setup before the C++ simulator is constructed).

``PHARE_H5_CHUNK_SIZE``
-------------------------

HDF5 dataset chunk size, in elements. Defaults to ``1024``. Affects diagnostic/restart HDF5 I/O.

``PHARE_SIM_MON``
-------------------

Enables periodic resource-usage monitoring during a simulation run. Truthy values: ``true``,
``1``, ``t`` (case-insensitive). Can be overridden per-run via ``Simulator.run(monitoring=...)``.

``PHARE_DRY_RUN``
-------------------

If ``1``, defaults ``Simulation.dry_run`` to ``True``, which makes ``Simulator.run()`` a no-op
(returns immediately without advancing). Used for CI dry-run validation of functional test scripts
without actually executing them.

``PYPHARE_LOG_LEVEL``
-----------------------

Sets the Python ``logging`` level used by PHARE's logger. One of ``INFO`` (default), ``WARNING``,
``ERROR``, ``DEBUG``.

``VIRTUAL_ENV``
-----------------

Standard Python venv marker, not PHARE-specific. If set, PHARE uses it to locate the venv's
``python3`` and pull its ``sys.path``, so an embedded/pybind11 Python (e.g. ``phare-exe``) can see
the venv's installed packages.

``OMPI_COMM_WORLD_RANK``, ``SLURM_PROCID``
---------------------------------------------

Not PHARE-specific - standard OpenMPI/Slurm launcher variables. Used only as a fallback
rank-detection heuristic (to decide whether a process should print) when MPI isn't initialized
yet.
