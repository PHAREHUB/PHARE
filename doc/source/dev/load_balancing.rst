.. _dev-load_balancing:

========================
Load balancing internals
========================

.. container:: breadcrumb

   You are here: :doc:`/dev/index` › Load balancing internals. User view:
   :doc:`/user_guide/load_balancing`.

**In one sentence:** PHARE fills a per-cell workload array on every
:term:`patch`; :term:`SAMRAI`'s partitioner uses it whenever it shares a
:term:`level` between ranks, and PHARE asks it to re-share level 0 between
:term:`regrids <regridding>`.

What the partitioner sees
-------------------------

.. plot::
   :include-source: False
   :caption: Illustrative numbers: the workload of each cell, summed per rank; SAMRAI moves the cut until the sums agree within tol.

   from phare_draw import *
   from matplotlib.patches import Rectangle
   counts = [5, 5, 5, 5, 10, 10, 20, 20]
   titles = ["Workload (mode nppc): particles per cell, level 0, 2 ranks",
             "After SAMRAI re-partitions level 0"]
   fig, axes = strip(2, titles, height=1.7)
   for ax, cut in zip(axes, [4, 6]):
       for i, c in enumerate(counts):
           col = C["L1"] if i < cut else C["L2"]
           ax.add_patch(Rectangle((i, 0), 1, 1, fc=col, alpha=0.18, ec=C["grey"], lw=0.8))
           text(ax, i + 0.5, 0.5, str(c))
       l0, l1 = sum(counts[:cut]), sum(counts[cut:])
       ax.add_patch(Rectangle((0, 0), cut, 1, fill=False, ec=C["L1"], lw=2.4))
       ax.add_patch(Rectangle((cut, 0), 8 - cut, 1, fill=False, ec=C["L2"], lw=2.4))
       text(ax, cut / 2, -0.4, f"rank 0: {l0}", color=C["L1"], bold=True)
       text(ax, cut + (8 - cut) / 2, -0.4, f"rank 1: {l1}", color=C["L2"], bold=True)
   text(axes[0], 4, -0.95, "max/avg − 1 = 60/40 − 1 = 0.5  >  tol = 0.05   →  re-partition", color=C["bad"])
   text(axes[1], 4, -0.95, "max/avg − 1 = 40/40 − 1 = 0  ≤  tol", color=C["ok"])
   for ax in axes:
       ax.plot([0, 8], [-1.25, 1.1], alpha=0)
   finish(fig, axes)

* ``LoadBalancerManager`` registers a cell variable with no ghosts and makes it the partitioner's
  workload (:src:`src/amr/load_balancing/load_balancer_manager.hpp#L23-L45`); without it SAMRAI counts cells.
* Hybrid levels: particles per cell (``nppc``,
  :src:`src/amr/load_balancing/concrete_load_balancer_strategy_nppc.hpp#L77-L81`) or 1 (``homogeneous``);
  MHD levels: always 1 (:src:`src/simulator/simulator.hpp#L275-L295`).
* Refreshed after each ``advanceLevel``, each new level, each restart load
  (:src:`src/amr/multiphysics_integrator.hpp#L338-L341`, :src:`src/amr/multiphysics_integrator.hpp#L558`).
* One ``CascadePartitioner``, tolerance ``tol``, used by SAMRAI's ``GriddingAlgorithm`` for every level
  (:src:`src/simulator/simulator.hpp#L257-L263`, :src:`src/amr/wrappers/integrator.hpp#L184-L186`).
  A new finest level is shared by cell count (SAMRAI uses the workload only for levels that already exist).

When level 0 is re-partitioned
------------------------------

.. plot::
   :include-source: False
   :caption: With auto, a balanced check doubles the wait, starting one check later; a rebalance resets it to 200.

   from phare_draw import *
   fig, (ax,) = strip(1, ["Level 0 checks with auto (next_rebalance=200, multiplier 2, max 1000)"], height=1.9)
   checks = [(0, False), (200, False), (600, False), (1400, True), (1600, False), (1800, False), (2200, False)]
   s = 1 / 200.0
   ax.plot([-0.3, 12.0], [0, 0], color=C["grey"], lw=1.2)
   for k in range(0, 2401, 400):
       ax.plot([k * s, k * s], [-0.08, 0.08], color=C["grey"], lw=1)
       text(ax, k * s, -0.38, str(k), fs=11)
   text(ax, 12.4, 0, "step", fs=11, ha="left")
   for (a, _), (b, _) in zip(checks, checks[1:]):
       arrow(ax, (a * s, 0.45), (b * s, 0.45), color=C["grey"], lw=1.2, ms=8, text=str(b - a), style="<|-|>")
   for st, reb in checks:
       mark(ax, "primal", st * s, 0, color=C["bad"] if reb else C["ok"], size=12)
   text(ax, 1400 * s, -0.85, "imbalanced: rebalance, wait back to 200", color=C["bad"], fs=11)
   ax.plot([0, 13], [-1.1, 0.8], alpha=0)
   finish(fig, [ax])

* ``Integrator::advance`` passes ``_should_rebalance_now()`` to ``advanceHierarchy(dt, …)``
  (:src:`src/amr/wrappers/integrator.hpp#L57-L64`): false unless tagging, hybrid level 0, and active
  (:src:`src/amr/wrappers/integrator.hpp#L204-L211`).
* ``every=N``: re-partition every ``N`` steps and at step 0 if ``on_init``, no measurement
  (``cadence_rebalance_check``, :src:`src/amr/wrappers/integrator.hpp#L106-L111`).
* ``auto``: ``CoarsestRebalanceCheck`` gathers the level 0 loads and tests ``max/avg − 1 > tol``
  (:src:`src/amr/wrappers/integrator.hpp#L113-L119`,
  :src:`src/amr/load_balancing/coarsest_rebalance_check.hpp#L35-L70`); chosen at
  :src:`src/amr/wrappers/integrator.hpp#L152-L153`.

The tricks
----------

One tolerance, two jobs
   **Problem.** The old level 0 test (``1 − min/max``) disagreed with SAMRAI's: level 0 was re-split forever.

   **Fix.** Same definition as SAMRAI, ``max/avg − 1 > tol`` (:src:`src/amr/load_balancing/coarsest_rebalance_check.hpp#L33-L47`).
   `#1308 <https://github.com/PHAREHUB/PHARE/pull/1308>`_ (load balancing cleanup).

Every rank agrees when to gather
   **Problem.** The load gather is collective: a rank that skips it blocks the others.

   **Fix.** The check counters are the same on every rank, so all decide together
   (:src:`src/amr/load_balancing/coarsest_rebalance_check.hpp#L16-L17`).

A private SAMRAI sum, redone
   **Problem.** ``CascadePartitioner::computeNonUniformWorkLoad`` is private.

   **Fix.** ``computeNonUniformWorkLoadForLevel0`` redoes the sum with the public
   ``BalanceUtilities`` (:src:`src/amr/wrappers/integrator.hpp#L90-L102`).

Not with refinement boxes
   **Problem.** Re-partitioning with user refinement boxes made SAMRAI abort in ``RefineSchedule``.

   **Fix.** ``Simulator::finalize_init`` refuses the combination (:src:`src/simulator/simulator.hpp#L264-L273`).

After a level 0 re-split, every level is stale
   **Problem.** SAMRAI reports levels (0, 0) changed, but finer levels read level 0.

   **Fix.** Every messenger is notified (:src:`src/amr/multiphysics_integrator.hpp#L372-L377`), see
   :doc:`/dev/messengers`.

Go deeper
---------

* :doc:`/dev/messengers` — what is rebuilt after a hierarchy change.
* :doc:`/dev/tagging_regridding` — when SAMRAI rebuilds levels.

Related
-------

* :doc:`/user_guide/load_balancing` — options and defaults.
* :doc:`/hpc/running_at_scale` — patch sizes and rank counts.
* :doc:`/dev/hpc_internals`
