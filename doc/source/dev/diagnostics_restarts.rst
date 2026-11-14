.. _dev-diagnostics_restarts:

================================
Diagnostics and restarts inside
================================

.. container:: breadcrumb

   You are here: :doc:`/dev/index` › Diagnostics and restarts. User view:
   :doc:`/user_guide/diagnostics`, :doc:`/user_guide/restarts`.

**In one sentence:** after each step, Python asks C++ to write a
:term:`restart` (one file per :term:`rank`, through :term:`SAMRAI`) and the
:term:`diagnostics` (HDF5 files shared by all ranks) when a requested time is
reached.

When a dump happens
-------------------

.. plot::
   :include-source: False
   :caption: Look at panel 2: in double precision the dump would fire at t = 0.2, one step before the requested 0.3.

   from phare_draw import *
   titles = ["The rule: dump at the step t where t_req − t < dt",
             "Double precision: 0.3 − 0.2 = 0.09999999999999998 < 0.1",
             "Cast to float first: 0.1 < 0.1 is false"]
   fig, axes = strip(3, titles, height=1.5)
   def axis(ax, dump_at, color, note):
       ax.plot([-0.3, 5.3], [0, 0], color=C["grey"], lw=1.2)
       for k in range(6):
           ax.plot([k, k], [-0.12, 0.12], color=C["ink"], lw=1.2)
           text(ax, k, -0.4, f"t = {k/10:.1f}" if k == 0 else f"{k/10:.1f}")
       ax.plot([3], [0.45], marker="v", ms=11, color=C["E"]); text(ax, 3.75, 0.5, "t_req = 0.3", color=C["E"], ha="left")
       mark(ax, "primal", dump_at, 0, color=color, size=14)
       text(ax, dump_at - 0.2, 0.45, note, color=color, bold=True, ha="right")
       ax.plot([-0.3, 5.3], [-0.65, 0.75], alpha=0)
   axis(axes[0], 3, C["ok"], "dump")
   axis(axes[1], 2, C["bad"], "dump, one step early")
   axis(axes[2], 3, C["ok"], "dump")
   finish(fig, axes)

* ``Simulator.initialize`` dumps once, ``Simulator.advance`` after every
  step; restarts first, then diagnostics
  (:src:`pyphare/pyphare/simulator/simulator.py#L130-L165`,
  :src:`pyphare/pyphare/simulator/simulator.py#L230-L238`).
* A requested time ``t_req`` is due when ``t_req − t < dt``; all passed
  times are consumed at once (``needsAction_``,
  :src:`src/diagnostic/diagnostic_manager.hpp#L124-L144`; restarts:
  :src:`src/restarts/restarts_manager.hpp#L74-L94`).
* Data are labelled with the step time ``t``, not ``t_req``.
* ``t`` is ``start + dt × step`` (constant ``dt``) or a compensated sum
  (adaptive ``dt``) (:src:`src/simulator/simulator.hpp#L577`,
  :src:`src/core/utilities/timestamps.hpp#L30-L100`): a dump is never
  missed or doubled, but could fire one step early (panel 2).

Writing in parallel (``phareh5``)
---------------------------------

.. plot::
   :include-source: False
   :caption: Rank 0 owns patches p0 and p1, rank 1 owns p2. Both create all three datasets, then each writes its own.

   from phare_draw import *
   from matplotlib.patches import Rectangle
   titles = ["Pass 1: every rank creates every dataset", "Pass 2: every rank writes only its own patches"]
   fig, axes = strip(2, titles, height=2.6)
   for ax, own in zip(axes, [False, True]):
       for r, (y, pats) in enumerate([(2.0, ["p0", "p1"]), (0.0, ["p2"])]):
           ax.add_patch(Rectangle((0, y - 0.35), 1.8, 0.7, fc=C["patch"], ec=C["patch_edge"], lw=1.5))
           text(ax, 0.9, y + 0.12, f"rank {r}", bold=True); text(ax, 0.9, y - 0.18, " ".join(pats), fs=11)
       for k, (name, y, owner) in enumerate([("p0", 2.0, 0), ("p1", 1.0, 0), ("p2", 0.0, 1)]):
           ax.add_patch(Rectangle((4.2, y - 0.3), 1.6, 0.6, fc=C["faint"], ec=C["ink"], lw=1.2))
           text(ax, 5.0, y, f"dataset {name}")
           for r, ry in enumerate([2.0, 0.0]):
               if own and r != owner:
                   continue
               arrow(ax, (1.85, ry), (4.15, y), color=C["L1"] if r == 0 else C["L2"], lw=1.6, ms=10)
       text(ax, 5.0, 2.65, "HDF5 file", bold=True)
   finish(fig, axes)

* HDF5 (with MPI) requires every rank to create every dataset, with the
  same arguments.
* Pass 1 gathers each dataset's path and size from all ranks
  (``create_data_set_per_mpi``, :src:`src/hdf5/detail/h5/h5_file.hpp#L138-L158`;
  attributes: :src:`src/hdf5/detail/h5/h5_file.hpp#L204-L258`); pass 2
  walks the patches and writes (:src:`src/diagnostic/detail/h5writer.hpp#L295-L364`).
* The writer is picked by ``format``: ``"phareh5"`` (default) or
  ``"pharevtkhdf"`` (:src:`src/diagnostic/diagnostics.hpp#L47-L66`).
  VTKHDF writes no ``particle``, ``meta``, ``info`` diagnostics, silently in
  a release build (:src:`src/diagnostic/detail/vtkh5_writer.hpp#L70-L86`).

Restarts
--------

.. code-block:: text

   <dir>/00000.30000/restore.000000/nodes.0000002/proc.0000000
         └ time,      └ always 0     └ number    └ rank
           11 chars, 5 decimals        of ranks

* Written by SAMRAI's ``RestartManager``, one file per rank; PHARE adds the
  patch data ids and the serialized simulation
  (:src:`src/restarts/detail/h5writer.hpp#L42-L60`).
* Path built in ``getRestartFileFullPath``
  (:src:`src/amr/wrappers/hierarchy.hpp#L75-L85`); on load, Python and
  SAMRAI open ``nodes.<current number of ranks>``: restart with the same
  number of ranks.
* Load: ``restart_time="auto"`` picks the latest saved time
  (:src:`pyphare/pyphare/pharein/restarts.py#L57-L74`); only ``cells``,
  ``dl``, ``max_nbr_levels`` are compared with the saved run
  (:src:`pyphare/pyphare/pharein/restarts.py#L167-L187`).
* ``keep_last`` deletes older time folders after a new restart
  (:src:`pyphare/pyphare/pharein/restarts.py#L17-L52`).

The tricks
----------

A rank with no patch still takes part
   **Problem.** A rank skipping an empty level would leave the others waiting in the collective calls.

   **Fix.** Every rank gets an entry, possibly empty, for every level up to the finest of all ranks
   (:src:`src/diagnostic/detail/h5writer.hpp#L320-L324`).

Wall-clock times come from rank 0
   **Problem.** ``elapsed_timestamps`` compared to each node's clock: ranks could disagree.

   **Fix.** ``mpi::unix_timestamp_now`` broadcasts rank 0's clock (:src:`src/mpi/mpi_utils.cpp#L93-L96`).

Open the restart file before the hierarchy exists
   **Problem.** SAMRAI's ``PatchHierarchy`` constructor checks whether the run is a restart.

   **Fix.** ``amr::Hierarchy`` inherits first from ``HierarchyRestarter``, which opens the file
   (:src:`src/amr/wrappers/hierarchy.hpp#L36-L73`, :src:`src/amr/wrappers/hierarchy.hpp#L112`).

Fine-level and emergency dumps
   **Problem.** Normal dumps happen only at level 0 step ends; a crash loses the state that caused it.

   **Fix.** ``fine_dump_lvl_max``: a refined level dumps all its diagnostics at its own sub-step ends,
   ignoring write timestamps (:src:`src/amr/multiphysics_integrator.hpp#L464-L480`,
   :src:`src/diagnostic/diagnostic_manager.hpp#L222-L231`). ``allow_emergency_dumps``: some errors
   dump all levels first (:src:`src/simulator/simulator.hpp#L607-L618`), only in builds where errors
   are exceptions (``devMode``, :doc:`/dev/hpc_internals`).

.. note:: Known issues

   * 3 or more MPI ranks: two identical runs can differ in the last bit; restarts inherit it:
     `#1526 <https://github.com/PHAREHUB/PHARE/issues/1526>`_ (MPI runs not bitwise).
   * Restarts write redundant data: `#664 <https://github.com/PHAREHUB/PHARE/issues/664>`_
     (restarts write redundant information).

Go deeper
---------

* :doc:`/dev/hpc_internals` — what is collective, errors across ranks.

Related
-------

* :doc:`/user_guide/diagnostics`, :doc:`/user_guide/restarts`
* :doc:`/dev/resources_manager` — patch data ids.
* :doc:`/hpc/running_at_scale` — I/O at scale.
