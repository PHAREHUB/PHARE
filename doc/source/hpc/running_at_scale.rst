================
Running at scale
================

**Goal:** run PHARE on many cores and avoid the usual traps: idle ranks,
parallel output, restarts.

.. contents:: On this page
   :local:
   :depth: 1


How PHARE runs in parallel
--------------------------

PHARE runs in parallel with :term:`MPI` only: you start one process (an MPI
:term:`rank`) per core, and there are no threads inside a rank. Each
:term:`level` of the mesh is cut into :term:`patches <patch>`. Each patch
belongs to one rank, which advances it and exchanges border values
(:term:`ghost cells <ghost cell>` and particles) with the ranks that own its neighbours.

.. raw:: html

   <svg viewBox="0 0 640 170" role="img" aria-label="Two mesh levels cut into patches, each patch coloured by the rank that owns it" style="max-width:100%;height:auto;font-family:sans-serif;font-size:14px">
     <text x="0" y="40" fill="currentColor">level 0</text>
     <rect x="80"  y="20" width="135" height="34" fill="#4c78a8" opacity="0.85"/>
     <rect x="215" y="20" width="135" height="34" fill="#f58518" opacity="0.85"/>
     <rect x="350" y="20" width="135" height="34" fill="#54a24b" opacity="0.85"/>
     <rect x="485" y="20" width="135" height="34" fill="#b279a2" opacity="0.85"/>
     <text x="0" y="100" fill="currentColor">level 1</text>
     <rect x="250" y="80" width="90"  height="34" fill="#4c78a8" opacity="0.85"/>
     <rect x="340" y="80" width="60"  height="34" fill="#f58518" opacity="0.85"/>
     <rect x="400" y="80" width="60"  height="34" fill="#54a24b" opacity="0.85"/>
     <rect x="460" y="80" width="40"  height="34" fill="#b279a2" opacity="0.85"/>
     <g fill="none" stroke="currentColor" stroke-width="1">
       <rect x="80" y="20" width="540" height="34"/>
       <rect x="250" y="80" width="250" height="34"/>
     </g>
     <rect x="80"  y="140" width="14" height="14" fill="#4c78a8"/><text x="100" y="152" fill="currentColor">rank 0</text>
     <rect x="180" y="140" width="14" height="14" fill="#f58518"/><text x="200" y="152" fill="currentColor">rank 1</text>
     <rect x="280" y="140" width="14" height="14" fill="#54a24b"/><text x="300" y="152" fill="currentColor">rank 2</text>
     <rect x="380" y="140" width="14" height="14" fill="#b279a2"/><text x="400" y="152" fill="currentColor">rank 3</text>
   </svg>

*Each rank owns patches on every level where it has work. The refined
level (level 1) covers only part of the domain, and its patches are shared
between ranks too.*

Who owns which patch is decided by load balancing, and changes during the
run: see :doc:`/user_guide/load_balancing`.

Launching
---------

.. code-block:: bash

   mpirun -n 64 python3 my_run.py

The same script runs on 1 or many ranks; you do not change the input. See
:doc:`/user_guide/running` for the other ways to start a run.

Patches and ranks
-----------------

A rank with no patch on a level has nothing to do on that level. With too
many ranks for the size of the problem, ranks wait for each other. Two
options set the size of patches, in cells of the level, per direction:

``smallest_patch_size`` (default 3)
   Patches are not cut smaller than this when a level is shared between
   ranks. A refined patch can still end up smaller: PHARE keeps it small
   rather than growing it across the domain boundary or over another
   patch.

``largest_patch_size`` (default: no limit)
   Patches larger than this are cut.

.. code-block:: python

   sim = ph.Simulation(
       # ...
       smallest_patch_size=20,
       largest_patch_size=100,
   )

The same values apply to every level.

.. admonition:: In preparation
   :class: in-preparation

   Recommended number of cells and particles per rank, and patch sizes, from
   measured runs. See also :doc:`scaling_results`.

Output at scale
---------------

* **Parallel HDF5 is required.** Output files are HDF5 files. With more
  than one rank, all ranks write into the same files. If your HDF5 library is not built for MPI, the
  run stops with ``HDF5 NOT PARALLEL!``.
* **Two output formats.** ``"phareh5"`` (default) and ``"pharevtkhdf"``
  (VTKHDF, the HDF5-based format of VTK, readable by ParaView). The VTKHDF
  writer does one write per rank and per level at each dump;
  `#1505 <https://github.com/PHAREHUB/PHARE/pull/1505>`_ (one write per
  rank) made it much faster at large rank counts. Its HDF5 chunk size can be tuned with
  ``PHARE_VTKHDF_CHUNK_BYTES`` (see :doc:`/user_guide/environment_variables`).
  It does not write particle diagnostics.
* **Log files.** By default each rank writes its own log file under
  ``.log/`` (``PHARE_LOG=RANK_FILES``). On thousands of ranks, set
  ``PHARE_LOG=NULL`` to discard them, or ``PHARE_LOG=CLI`` to print to the
  terminal (see :doc:`/user_guide/environment_variables`).

Restarts at scale
-----------------

* **Restart with the same number of ranks.** A restart is saved as one file
  per rank, in a folder named after the number of ranks. Restarting with a
  different ``mpirun -n`` fails to find the files.
* **Save before the job's time limit.** ``elapsed_timestamps`` asks for a
  restart at given wall-clock times since the start (in seconds, or
  ``datetime.timedelta``), whatever the simulation time:

  .. code-block:: python

     from datetime import timedelta

     sim = ph.Simulation(
         # ...
         restart_options={"dir": "checkpoints", "mode": "overwrite",
                          "elapsed_timestamps": [timedelta(hours=11, minutes=30)]},
     )

  All ranks use the clock of rank 0, so they decide together.
* **Disk space.** ``keep_last=N`` in ``restart_options`` deletes older
  restart folders and keeps the ``N`` most recent.

See :doc:`/user_guide/restarts` for all restart options.

When something fails
--------------------

When a rank hits an error during initialization or a time step, what
happens depends on how PHARE was built. (Errors at setup, such as a wrong
option, and while writing output are ordinary Python exceptions in any
build.)

* **Default build:** the rank prints ``MPI ABORT: PHARE ERROR:`` followed
  by the place in the code and the message, and the whole MPI job is
  aborted at once. The message goes to that rank's log file under
  ``.log/`` (with the default ``PHARE_LOG=RANK_FILES``). With
  ``PHARE_LOG=NULL`` this message is lost: use ``CLI`` when you need it.
* **Built with** ``-DuseExceptionsInsteadOfMPIAbort=ON`` **(or**
  ``-DdevMode=ON`` **):** the failing rank prints its message on the
  standard error. At the end of the step all ranks check together whether
  any of them failed; they then close the diagnostic files and the script
  stops with exit code 1.

Common mistakes
---------------

* **Restarting on a different number of ranks** — the run stops at setup
  with ``Unable to open file .../restore.000000/nodes.0000001/proc.0000000``
  (here: saved on 2 ranks, restarted on 1) — restart with the same
  ``mpirun -n``.
* **Serial HDF5 on a cluster** — ``HDF5 NOT PARALLEL!`` — load or build an
  MPI-enabled HDF5 and rebuild PHARE against it.
* **Many more ranks than patches** — no speed-up, ranks idle — use fewer
  ranks or smaller ``largest_patch_size``.

.. seealso:: What happens inside during parallel output and load balancing:
   :doc:`/dev/hpc_internals`.
