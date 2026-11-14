.. _dev-hpc_internals:

=============
HPC internals
=============

.. container:: breadcrumb

   You are here: :doc:`/dev/index` › HPC internals. User view: :doc:`/hpc/running_at_scale`.

**In one sentence:** PHARE is :term:`MPI`-only (one process per core, no threads); this page lists
what every :term:`rank` must do together, how errors cross ranks, and where the timers sit.

What every rank must do together
--------------------------------

.. plot::
   :include-source: False
   :caption: Illustrative timings: at each dashed line every rank must arrive; the slowest rank sets the pace.

   from phare_draw import *
   from matplotlib.patches import Rectangle
   fig, (ax,) = strip(1, ["One time step on 3 ranks: every vertical line is a collective call"], height=3.2)
   # (name, x) of the collective points, in step order
   cols = [("dt\nmin", 1.0), ("level 0\ngather", 2.4), ("errors?\n(push)", 6.0), ("errors?\n(step end)", 8.4), ("dump", 10.0)]
   work = {0: [(1.0, 2.4, 0.8), (2.4, 6.0, 2.5), (6.0, 8.4, 1.6), (8.4, 10.0, 1.2)],
           1: [(1.0, 2.4, 0.9), (2.4, 6.0, 2.2), (6.0, 8.4, 1.9), (8.4, 10.0, 1.2)],
           2: [(1.0, 2.4, 1.3), (2.4, 6.0, 3.5), (6.0, 8.4, 2.3), (8.4, 10.0, 1.2)]}
   for r in range(3):
       y = 2 - r
       text(ax, -0.15, y, f"rank {r}", ha="right", bold=True)
       for x0, x1, busy in work[r]:
           ax.add_patch(Rectangle((x0 + 0.05, y - 0.22), busy, 0.44, fc=C["patch"], ec=C["patch_edge"], lw=1.2))
           if x0 + 0.05 + busy < x1 - 0.05:
               ax.add_patch(Rectangle((x0 + 0.05 + busy, y - 0.22), x1 - 0.1 - x0 - busy, 0.44,
                                      fc=C["faint"], ec=C["grey"], lw=0.8, hatch="///"))
       ax.add_patch(Rectangle((0.15, y - 0.22), 0.8, 0.44, fc=C["patch"], ec=C["patch_edge"], lw=1.2))
   for name, x in cols:
       ax.plot([x, x], [-0.45, 2.45], color=C["bad"], lw=1.8, ls="--")
       text(ax, x, 2.95, name, color=C["bad"], fs=11)
   text(ax, 4.2, -0.75, "advance + ghost exchanges", fs=11)
   text(ax, 9.2, -0.75, "Python", fs=11)
   ax.add_patch(Rectangle((0.2, -1.25), 0.4, 0.25, fc=C["patch"], ec=C["patch_edge"])); text(ax, 0.75, -1.12, "work", ha="left")
   ax.add_patch(Rectangle((2.2, -1.25), 0.4, 0.25, fc=C["faint"], ec=C["grey"], hatch="///")); text(ax, 2.75, -1.12, "waiting for the slowest rank", ha="left")
   ax.plot([-1.2, 10.4], [-1.4, 3.4], alpha=0)
   finish(fig, [ax])

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Collective call (step order)
     - Code
   * - adaptive ``dt``: ``mpi::min``
     - ``MultiPhysicsIntegrator::computeStableDt`` (:src:`src/amr/multiphysics_integrator.hpp#L435-L460`)
   * - level 0 load gather (``auto`` only, when due)
     - ``tolerance_rebalance_check`` (:src:`src/amr/wrappers/integrator.hpp#L113-L119`), :doc:`load_balancing`
   * - ``mpi::any_errors`` (exception builds)
     - after the particle push (:src:`src/amr/solvers/solver_ppc.hpp#L614-L615`), the MHD level advance
       (:src:`src/amr/solvers/solver_mhd.hpp#L381-L382`), end of ``initialize`` and of each step
       (:src:`src/simulator/simulator.hpp#L542`, :src:`src/simulator/simulator.hpp#L596`)
   * - wall-clock checks: rank 0's clock
     - ``mpi::unix_timestamp_now`` (:src:`src/mpi/mpi_utils.cpp#L93-L96`)
   * - restart dump: barrier; diagnostics dump: HDF5 metadata
     - :src:`src/restarts/detail/h5writer.hpp#L59`, :doc:`diagnostics_restarts`
   * - Python data access
     - ``DataWrangler::sync`` gathers every patch (:src:`src/python3/data_wrangler.hpp#L75-L102`):
       reading data on rank 0 only would hang.

PHARE's own MPI helpers are in ``src/mpi/mpi_utils.hpp``
(:src:`src/mpi/mpi_utils.hpp#L22-L47`, :src:`src/mpi/mpi_utils.hpp#L115-L120`): ``collect``, ``max``,
``min``, ``any``, ``barrier``, ``all_get_from_rank_0``.

Errors across ranks
-------------------

* **Default build:** ``PHARE_LOG_ERROR`` prints ``MPI ABORT: PHARE ERROR: …`` and calls ``MPI_Abort``
  (``mpi::log_error``, :src:`src/mpi/mpi_utils.cpp#L65-L71`); ``any_errors`` is always false.
* **Exception builds** (``useExceptionsInsteadOfMPIAbort``, forced by ``devMode``, :src:`res/cmake/def.cmake#L25-L31`): the
  failing rank prints to standard error; at the next ``any_errors`` all ranks close diagnostics and throw
  (:src:`src/simulator/simulator.hpp#L596-L602`); Python prints on rank 0 and exits with code 1 while
  ``exit_on_exception`` is true (the default, :src:`pyphare/pyphare/simulator/simulator.py#L314-L318`).

Timers
------

* ``PHARE_LOG_SCOPE(level, "Name::part")`` is kept when ``level ≤ PHARE_LOG_LEVEL``, a compile-time macro,
  default 0 (:src:`src/core/logger.hpp#L10-L16`).
* phlop backend (``-DwithPhlop=ON``): macros in :src:`src/core/utilities/logger/logger_defaults.hpp#L17-L39`;
  timing starts if ``PHARE_SCOPE_TIMING`` is ``1`` or ``true`` (:src:`src/amr/samrai.cpp#L25-L28`).
* Caliper backend (``-DwithCaliper=ON``): ``PHARE::scope_log`` (:src:`src/core/logger.hpp#L31-L60`).

The tricks
----------

Decide collectives from shared data only
   **Problem.** Ranks that disagree on entering a collective call wait forever.

   **Fix.** Decisions use data all ranks share: step counters, rank 0's clock, gathered values.

Check for errors inside the step
   **Problem.** In exception builds a failing rank leaves its step; the others would wait at the next exchange.

   **Fix.** ``any_errors`` right after the parts that can fail, then throw everywhere
   (:src:`src/amr/solvers/solver_ppc.hpp#L610-L615`).

Data layout and hot loops
-------------------------

.. admonition:: In preparation
   :class: in-preparation

   Memory layout of particles and fields, the hot loops (push, interpolation, deposition) and their share
   of a step. Until then: :doc:`particles_ions`, :doc:`fields`, :doc:`ppc_order`, :doc:`/hpc/memory`.

Go deeper
---------

* :doc:`load_balancing` · :doc:`diagnostics_restarts` · :doc:`messengers`

Related
-------

* :doc:`/hpc/measuring_performance`
* :doc:`working_on_phare` — build options.
