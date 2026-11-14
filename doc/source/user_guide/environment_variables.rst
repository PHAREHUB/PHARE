=====================
Environment variables
=====================

**Goal:** know which environment variables change how PHARE runs, and set them
without editing your script.

An environment variable is set in the shell, either for one run
(``PHARE_LOG=CLI python3 my_script.py``) or for the whole session
(``export PHARE_LOG=CLI``). PHARE reads the ones below, on both the Python and
the C++ side.

.. list-table::
   :header-rows: 1
   :widths: 24 56 20

   * - Variable
     - What it does
     - Default
   * - ``PHARE_LOG``
     - Where the C++ text output goes: ``RANK_FILES`` (one file per MPI rank in
       ``.log/``), ``DATETIME_FILES`` (one file per rank named with the start
       date and time), ``CLI`` (the terminal), ``NULL`` (discarded). Any other
       value is an error.
     - ``RANK_FILES`` when running through Python (``Simulator`` log_to_file on);
       no redirection if unset in a pure C++ run
   * - ``PHARE_DRY_RUN``
     - If ``1``, sets ``Simulation(dry_run=...)`` to ``True``: the script is
       built and checked, but ``Simulator.run()`` does nothing.
     - ``0`` (off)
   * - ``PHARE_SIM_MON``
     - Turns on resource monitoring during ``Simulator.run()``. Accepted
       (case-insensitive): ``true``, ``1``, ``t``. An explicit
       ``run(monitoring=...)`` takes precedence. Needs the ``phlop`` package.
     - off
   * - ``PHARE_SCOPE_TIMING``
     - Turns on per-function timing, written to ``.phare/timings/rank.<N>.txt``
       (one file per rank). Needs a build with ``phlop`` and the ``phlop``
       package. Use ``1`` or ``true``.
     - off
   * - ``PHARE_VTKHDF_CHUNK_BYTES``
     - Size in bytes of the HDF5 chunks used for field data in ``pharevtkhdf``
       diagnostics.
     - ``1048576`` (1 MiB)
   * - ``PYPHARE_LOG_LEVEL``
     - Level of the messages of PHARE's Python logger: ``DEBUG``, ``INFO``,
       ``WARNING`` or ``ERROR`` (capitals).
     - ``INFO``
   * - ``VIRTUAL_ENV``
     - Not set by you: the Python virtual environment sets it. PHARE uses it to
       find the packages of the environment (needed by the embedded Python of
       ``phare-exe``).
     - set by the venv
   * - ``OMPI_COMM_WORLD_RANK``, ``PMIX_RANK``, ``PMI_RANK``, ``SLURM_PROCID``
     - Not PHARE variables: they are set by ``mpirun`` / ``srun``. PHARE reads
       them only to know whether a process is rank 0 (so that only one process
       prints) before MPI is started.
     - set by the launcher

Examples
--------

See the log in the terminal instead of in ``.log/``:

.. code-block:: bash

   PHARE_LOG=CLI python3 my_script.py

Check a script without running it (see :doc:`running`):

.. code-block:: bash

   PHARE_DRY_RUN=1 python3 my_script.py

Monitor CPU and memory every run, without touching the script (the sampling
interval is then 1 second):

.. code-block:: bash

   PHARE_SIM_MON=1 python3 my_script.py

For a chosen interval, use ``Simulator.run(monitoring=10)`` in the script.

Fewer, larger chunks for big VTKHDF fields (4 MiB):

.. code-block:: bash

   PHARE_VTKHDF_CHUNK_BYTES=4194304 python3 my_script.py

Notes on some variables
-----------------------

``PHARE_VTKHDF_CHUNK_BYTES``
   A chunk holds whole rows (all components of a point). The other VTKHDF
   datasets (time steps, offsets and AMR boxes) are small and use fixed chunks
   of 1024 rows.

.. _phare_scope_timing:

``PHARE_SCOPE_TIMING``
   The Python side also accepts ``t`` and capital letters; the C++ side, which
   writes the files, only accepts ``1`` or ``true``: use one of those. The
   Python side also prepares the ``.phare/timings`` folder. See
   :doc:`/hpc/measuring_performance` for configuring and reading the timings.

``PHARE_LOG``
   Python creates the ``.log/`` folder on rank 0 unless you chose ``CLI`` or
   ``NULL``. See :doc:`running`.
