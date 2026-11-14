.. _dev-diagnostics_restarts:

================================
Diagnostics and restarts inside
================================

.. container:: breadcrumb

   You are here: :doc:`/dev/index` ›
   Diagnostics and restarts. User view: :doc:`/user_guide/diagnostics`,
   :doc:`/user_guide/restarts`.

**In one sentence:** after each time step, Python asks C++ to write a
:term:`restart` and :term:`diagnostics` if a requested time has been
reached; diagnostics go through one *writer* per output format that writes
HDF5 files shared by all :term:`ranks <rank>`, restarts go through
:term:`SAMRAI`'s restart manager, one file per rank.

The picture
-----------

.. mermaid::

   flowchart TD
     I["Simulator.initialize / advance<br/>(Python)"] --> D["Simulator.dump()"]
     D --> R["1. restarts: RestartsManager<br/>then SAMRAI RestartManager<br/>(one file per rank)"]
     D --> G["2. diagnostics:<br/>DiagnosticsManager"]
     G --> W{"format"}
     W -- "phareh5" --> H5["h5::H5Writer<br/>all diagnostic types"]
     W -- "pharevtkhdf" --> VT["vtkh5::H5Writer<br/>fields only"]

How it works
------------

#. **Who triggers a dump.** Python's ``Simulator.initialize`` dumps once
   before the first step, and ``Simulator.advance`` dumps after each step
   (:src:`pyphare/pyphare/simulator/simulator.py#L130-L165`). ``dump()``
   does restarts first, then diagnostics
   (:src:`pyphare/pyphare/simulator/simulator.py#L230-L238`).
#. **Is it time?** ``DiagnosticsManager`` keeps, per diagnostic, the index
   of the next requested time; ``RestartsManager`` keeps one for its single
   restart block. At current time ``t``
   with step ``dt``, a requested time ``t_req`` is due when
   ``t_req - t < dt`` (``needsAction_``,
   :src:`src/diagnostic/diagnostic_manager.hpp#L124-L144`). All requested
   times already passed are consumed at once. Data are labelled with the
   actual time ``t``, not ``t_req``.
#. **Which diagnostics exist.** Hybrid runs know ``fluid``, ``electromag``,
   ``particle``, ``meta``, ``info``; MHD runs know ``mhd``, ``meta``,
   ``electromag``. Each user block is named type + number (``fluid0``,
   ``fluid1``...) (``registerDiagnostics``,
   :src:`src/diagnostic/diagnostic_manager.hpp#L28-L60`).
#. **The writer.** ``DiagnosticsManagerResolver`` picks the writer from
   ``format``: ``"phareh5"`` (default) or ``"pharevtkhdf"``
   (:src:`src/diagnostic/diagnostics.hpp#L47-L66`). The writer walks the
   hierarchy patch by patch and calls one *type writer* per diagnostic type
   (``H5Writer::dump``, :src:`src/diagnostic/detail/h5writer.hpp#L245-L273`).
   The VTKHDF writer (VTKHDF is the HDF5-based file layout of VTK, read by
   ParaView) maps ``particle``, ``meta`` and ``info`` to a writer that does
   nothing, so these diagnostics produce no output in that format, with no
   warning in a release build
   (:src:`src/diagnostic/detail/vtkh5_writer.hpp#L70-L86`).
#. **Restarts.** ``RestartsManager::dump`` asks the restart writer to call
   SAMRAI's ``RestartManager::writeRestartFile``, then adds the patch data
   ids and the serialized simulation to the rank's file
   (``restarts::h5::Writer::dump``,
   :src:`src/restarts/detail/h5writer.hpp#L42-L60`). Python then deletes
   older restart directories if ``keep_last`` is set
   (:src:`pyphare/pyphare/pharein/restarts.py#L17-L52`).
#. **Loading a restart.** Python finds the directory for ``restart_time``
   (``"auto"`` = latest), checks that ``cells``, ``dl`` and
   ``max_nbr_levels`` match the saved run
   (:src:`pyphare/pyphare/pharein/restarts.py#L167-L187`), and passes the
   saved patch data ids and the path to C++
   (:src:`pyphare/pyphare/pharein/initialize/general.py#L169-L192`).

The tricks
----------

Every rank creates every dataset
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** With parallel HDF5 (the HDF5 file library built on MPI-IO),
creating a dataset or attribute changes the
file's metadata, and HDF5 requires every rank to do it, with the same
arguments. But each rank only knows its own patches.

**Fix.** Writing is done in two passes. First each rank collects the path
and size of each dataset it will write; ``create_data_set_per_mpi`` gathers
all paths and sizes from all ranks and every rank creates all datasets
(:src:`src/hdf5/detail/h5/h5_file.hpp#L138-L158`; attributes:
``write_attributes_per_mpi``, :src:`src/hdf5/detail/h5/h5_file.hpp#L204-L240`).
Then each rank writes the data of its own patches
(``initializeDatasets_`` and ``writeDatasets_``,
:src:`src/diagnostic/detail/h5writer.hpp#L295-L364`).

.. mermaid::

   sequenceDiagram
     participant R0 as rank 0
     participant R1 as rank 1
     participant F as HDF5 file
     Note over R0,R1: pass 1: list own datasets
     R0->>R1: gather paths, sizes
     R1->>R0: gather paths, sizes
     R0->>F: create ALL datasets
     R1->>F: create ALL datasets
     Note over R0,R1: pass 2: write own patches
     R0->>F: write own data
     R1->>F: write own data

A rank with no patch still takes part
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** A rank may own no patch on some level. If it skipped that level,
the other ranks would wait forever in the collective calls.

**Fix.** The writer computes the finest level over all ranks
(``mpi::max``) and gives every rank an entry, possibly empty, for every
level up to it (:src:`src/diagnostic/detail/h5writer.hpp#L320-L324`,
:src:`src/diagnostic/detail/h5writer.hpp#L355-L359`).

Wall-clock times come from rank 0
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** Diagnostics and restarts can be requested at elapsed wall-clock
times (``elapsed_timestamps``). Each node has its own clock: ranks could
disagree on whether the time has come, and only some would enter the
collective write.

**Fix.** ``mpi::unix_timestamp_now`` reads rank 0's clock and broadcasts it
(:src:`src/mpi/mpi_utils.cpp#L93-L96`), so all ranks decide together.

Time comparison in single precision
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** Times like 0.1 are not exact in binary. One step before a
requested time, ``t_req - t`` should equal ``dt`` exactly, but it can come
out a hair below ``dt``, and the dump would happen one step early. (A dump
cannot be missed or doubled: the index of the next requested time always
moves forward, and ``t`` itself is computed as ``start + dt * step`` for a
constant step, with compensated summation for an adaptive one,
:src:`src/core/utilities/timestamps.hpp#L30-L66`.)

**Fix.** ``needsAction_`` compares ``t_req - t`` and ``dt`` after casting
both to ``float``, which drops those last digits, so the two compare equal
(:src:`src/diagnostic/diagnostic_manager.hpp#L124-L130`; same in
:src:`src/restarts/restarts_manager.hpp#L74-L80`).

Open the restart file before the hierarchy exists
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** SAMRAI's ``PatchHierarchy`` constructor checks whether the run
starts from a restart, so the restart file must already be open.

**Fix.** ``amr::Hierarchy`` inherits first from ``HierarchyRestarter``, whose
constructor registers the saved patch data ids and opens the file. Base
classes are built first, so this happens before the ``PatchHierarchy`` base
(:src:`src/amr/wrappers/hierarchy.hpp#L36-L73`,
:src:`src/amr/wrappers/hierarchy.hpp#L112`).

Restarts need the same number of ranks
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem / constraint.** SAMRAI writes one file per rank, at
``<dir>/<time>/restore.000000/nodes.<number of ranks>/proc.<rank>``
(``getRestartFileFullPath``,
:src:`src/amr/wrappers/hierarchy.hpp#L75-L85`). On load, both Python
(reading the saved simulation through ``getRestartFileFullPath``,
:src:`src/python3/cpp_etc.cpp#L83`) and SAMRAI
(:src:`src/amr/wrappers/hierarchy.hpp#L58-L59`) open the file for the
*current* number of ranks. A run restarted with a different number of ranks
stops at setup with "Unable to open file
.../restore.000000/nodes.0000004/proc.0000002" (here: saved on another
number of ranks, reloaded on 4).

``<time>`` is written with 11 characters and 5 decimals, e.g.
``00001.50000`` (``restartFilePathForTime``,
:src:`src/amr/wrappers/hierarchy.hpp#L122-L128`).

Dumps on fine levels and emergency dumps
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** Finer levels take several sub-steps per coarse step; normal
dumps only happen at coarse step ends. And when a run crashes, the state
that caused it is lost.

**Fix.** With ``fine_dump_lvl_max``, a ``pre_advance`` functor dumps a single
level at its own sub-step ends (``MultiPhysicsIntegrator::dump_``,
:src:`src/amr/multiphysics_integrator.hpp#L464-L480`; set up in
``Simulator::diagnostics_init``, :src:`src/simulator/simulator.hpp#L230-L253`).
With ``allow_emergency_dumps``, some known errors trigger a dump of all
levels before the exception is raised (``handle_dictionary_exception``,
:src:`src/simulator/simulator.hpp#L607-L618`). These errors are thrown
only after all ranks agree that one of them failed (``mpi::any_errors``,
e.g. :src:`src/amr/solvers/solver_ppc.hpp#L610-L615`). On any error, the
diagnostics manager is destroyed first so HDF5 files are flushed and closed
(:src:`src/simulator/simulator.hpp#L596-L602`).

Both emergency dumps and this clean stop only exist when errors are
exceptions, i.e. builds with
``useExceptionsInsteadOfMPIAbort`` (forced on by ``devMode``,
:src:`res/cmake/def.cmake#L25-L31`). In a default build,
``PHARE_LOG_ERROR`` calls ``MPI_Abort`` right away
(``mpi::log_error``, :src:`src/mpi/mpi_utils.cpp#L65-L71`) and
``mpi::any_errors`` always returns false
(:src:`src/mpi/mpi_utils.cpp#L55-L61`).

.. note:: Known issues

   * With 3 or more MPI ranks, two identical runs can differ in the last
     bit, and restarts inherit this: `#1526 <https://github.com/PHAREHUB/PHARE/issues/1526>`_ (MPI runs not bitwise).
   * Restarts write some redundant data: `#664 <https://github.com/PHAREHUB/PHARE/issues/664>`_ (restarts write redundant
     information).

Go deeper
---------

* :doc:`/dev/hpc_internals` — what is collective and what it costs.

Related
-------

* :doc:`/user_guide/diagnostics`, :doc:`/user_guide/restarts`
* :doc:`/dev/resources_manager` — patch data ids.
* :doc:`/hpc/running_at_scale` — I/O at scale.
