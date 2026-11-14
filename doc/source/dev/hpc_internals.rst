.. _dev-hpc_internals:

=============
HPC internals
=============

.. container:: breadcrumb

   You are here: :doc:`/dev/index` ›
   HPC internals. User view: :doc:`/hpc/running_at_scale`.

**In one sentence:** PHARE is :term:`MPI`-only (one process per core, no
threads); this page lists what every :term:`rank` must do together, how errors cross ranks,
and where the timers sit.

The picture
-----------

.. mermaid::

   flowchart TD
     DT["adaptive dt?<br/>mpi::min of stable steps"] --> LB["level 0 check due?<br/>gather loads"]
     LB --> ADV["advance my patches,<br/>exchange ghosts (SAMRAI schedules);<br/>any_errors after each push"]
     ADV --> ERR["mpi::any_errors()<br/>at end of step"]
     ERR --> DUMP["dump due? (Python)<br/>collective HDF5 writes"]
     click ADV href "messengers.html" "Messengers and schedules"
     click LB href "load_balancing.html" "Load balancing internals"
     click DUMP href "diagnostics_restarts.html" "Diagnostics and restarts"

How it works
------------

#. **One rank per core.** There is no thread parallelism in ``src/``: all
   parallelism is MPI, through SAMRAI (patches distributed over ranks,
   ghost exchanges by schedules) and a few PHARE helpers in
   ``src/mpi/mpi_utils.hpp`` (:src:`src/mpi/mpi_utils.hpp#L22-L47`,
   :src:`src/mpi/mpi_utils.hpp#L115-L120`).
#. **PHARE's own MPI helpers.** ``mpi::collect`` (gather a value from every
   rank to every rank), ``mpi::max``, ``mpi::min``, ``mpi::any``,
   ``mpi::barrier``, ``mpi::all_get_from_rank_0`` (compute on rank 0, give
   the result to all).
#. **Timers.** Code is instrumented with
   ``PHARE_LOG_SCOPE(level, "Name::part")``. A build keeps the scopes whose
   level is at most ``PHARE_LOG_LEVEL``, a compile-time macro (default 0 =
   none: add ``-DPHARE_LOG_LEVEL=1`` (to 3) to ``CMAKE_CXX_FLAGS``)
   (:src:`src/core/logger.hpp#L10-L16`,
   :src:`src/core/utilities/logger/logger_defaults.hpp#L17-L39`). The
   backend is phlop (``-DwithPhlop=ON``) or the Caliper profiling library
   (``-DwithCaliper=ON``). With phlop, timing starts only if
   ``PHARE_SCOPE_TIMING`` is ``1`` or ``true`` (:src:`src/amr/samrai.cpp#L25-L28`).

What every rank must call together
-----------------------------------

A *collective* call must be made by all ranks, in the same order; if one
rank skips it, the others wait forever. These are the main ones PHARE
adds on top of SAMRAI's own, in the order they occur in a step:

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Where
     - Why it is collective
   * - ``MultiPhysicsIntegrator::computeStableDt`` (adaptive time step,
       before the step)
     - ``mpi::min`` of the local stable steps
       (:src:`src/amr/multiphysics_integrator.hpp#L435-L460`).
   * - Level 0 load check (start of ``Integrator::advance``)
     - ``mpi::collect`` of each rank's load, when a check is due
       (:src:`src/amr/wrappers/integrator.hpp#L57-L64`, :doc:`load_balancing`).
   * - ``mpi::any_errors`` (exception builds)
     - all-reduce of the error flags: after each particle push
       (:src:`src/amr/solvers/solver_ppc.hpp#L614-L615`), after each MHD
       level advance (:src:`src/amr/solvers/solver_mhd.hpp#L381-L382`), and
       at the end of the step (:src:`src/simulator/simulator.hpp#L596`).
   * - Elapsed-time checks (dump decision)
     - ``mpi::unix_timestamp_now`` broadcasts rank 0's clock
       (:src:`src/mpi/mpi_utils.cpp#L93-L96`).
   * - Restart dump
     - ends with ``mpi::barrier``
       (:src:`src/restarts/detail/h5writer.hpp#L59`).
   * - Diagnostics dump (``H5Writer``)
     - HDF5 metadata: every rank creates every dataset
       (:doc:`diagnostics_restarts`).
   * - Python data access (``DataWrangler::sync``)
     - ``mpi::collect`` of every patch to every rank
       (:src:`src/python3/data_wrangler.hpp#L75-L102`): a ``post_advance``
       that reads data on rank 0 only would hang.

The tricks
----------

Decide collectives from shared data only
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** If ranks decide differently whether to enter a collective call
(a dump, a load gather), some enter and wait forever.

**Fix.** Every such decision is computed from data all ranks share: step
counters (``CoarsestRebalanceCheck``), rank 0's clock
(``mpi::unix_timestamp_now``), gathered values. Never from local data
alone.

Check for errors inside the step
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** In exception builds, a rank that fails leaves its step early;
the others would go on to the next ghost exchange and wait for it.

**Fix.** The solvers call ``mpi::any_errors()`` right after the parts that
can fail (particle push, MHD level advance) and throw on every rank at
once (:src:`src/amr/solvers/solver_ppc.hpp#L610-L615`).

Errors across ranks
-------------------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Build
     - What happens on an error
   * - default
     - ``PHARE_LOG_ERROR`` prints ``MPI ABORT: PHARE ERROR: <file:line>
       <message>`` and calls ``MPI_Abort``: the whole job stops at once
       (``mpi::log_error``, :src:`src/mpi/mpi_utils.cpp#L65-L71`).
   * - ``useExceptionsInsteadOfMPIAbort`` (forced by ``devMode``)
     - the failing rank records the error and prints it on standard error;
       at the next ``mpi::any_errors`` all ranks see it, close diagnostics
       and throw (:src:`src/simulator/simulator.hpp#L596-L602`,
       :src:`res/cmake/def.cmake#L25-L31`); Python prints the message on
       rank 0 and the script exits with code 1
       (:src:`pyphare/pyphare/simulator/simulator.py#L314-L318`).

Data layout and hot loops
-------------------------

.. admonition:: In preparation
   :class: in-preparation

   Memory layout of particles and fields, the hot loops of a hybrid step
   (push, interpolation, deposition) and their measured share of a step.
   Until then see :doc:`particles_ions`, :doc:`fields` and
   :doc:`ppc_order`.

Where memory goes
-----------------

.. admonition:: In preparation
   :class: in-preparation

   Bytes per particle and per cell, ghost overheads, SAMRAI metadata;
   user-side view in :doc:`/hpc/memory`.

Go deeper
---------

* :doc:`load_balancing`
* :doc:`diagnostics_restarts`
* :doc:`messengers`

Related
-------

* :doc:`/hpc/measuring_performance`
* :doc:`working_on_phare` — build options.
