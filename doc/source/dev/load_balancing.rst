.. _dev-load_balancing:

========================
Load balancing internals
========================

.. container:: breadcrumb

   You are here: :doc:`/dev/index` ›
   Load balancing. User view: :doc:`/user_guide/load_balancing`.

**In one sentence:** PHARE fills a per-cell "workload" array on every
:term:`patch` and hands it to :term:`SAMRAI`'s partitioner, which uses it each
time a :term:`level` is rebuilt; on top of that, PHARE asks SAMRAI to
re-partition level 0 between :term:`regrids <regridding>` when its load has
drifted (``auto``) or every N steps (``every``).

The picture
-----------

**Set-up, once** (which load SAMRAI uses when it shares patches):

.. mermaid::

   flowchart TD
     PY["pharein.LoadBalancer"] --> DET["LoadBalancerDetails"]
     DET --> CP["SAMRAI CascadePartitioner<br/>tolerance = tol"]
     DET --> LBM["LoadBalancerManager<br/>per-cell workload array"]
     LBM -- "workload index" --> CP
     LBM --> EST["one estimator per level<br/>nppc or homogeneous"]
     CP --> GA["GriddingAlgorithm<br/>every regrid, every level"]
     click PY href "../user_guide/load_balancing.html" "User options"
     click GA href "tagging_regridding.html" "Tagging and regridding"

**Every time step** (the extra level 0 check):

.. mermaid::

   flowchart TD
     ADV["Integrator::advance"] --> CHK["_should_rebalance_now<br/>(every N steps, or<br/>CoarsestRebalanceCheck)"]
     CHK --> AH["advanceHierarchy(dt, rebalance)"]
     AH -- "if rebalance" --> RP["SAMRAI re-partitions level 0"]
     RP --> RHC["resetHierarchyConfiguration:<br/>every messenger notified"]
     click RHC href "messengers.html" "Messengers"

How it works
------------

#. **Options.** ``pharein`` writes the options under
   ``simulation/AMR/loadbalancing``. ``active`` is true when the user
   declared a ``LoadBalancer`` or uses tagging
   (``populate_load_balancer``, :src:`pyphare/pyphare/pharein/initialize/general.py#L83-L102`).
   C++ reads them into ``LoadBalancerDetails``
   (:src:`src/amr/load_balancing/load_balancer_details.hpp#L37-L51`).
#. **One partitioner.** ``Simulator::finalize_init`` builds one SAMRAI
   ``CascadePartitioner`` with ``flexible_load_tolerance = tol``
   (:src:`src/simulator/simulator.hpp#L257-L263`) and passes it to SAMRAI's
   ``GriddingAlgorithm`` (:src:`src/amr/wrappers/integrator.hpp#L184-L186`).
   Every time SAMRAI builds or rebuilds a level, this partitioner decides
   which rank owns which patch.
#. **A workload array.** ``LoadBalancerManager`` registers a cell-centred
   SAMRAI variable, ``LoadBalancerVariable``, with no ghost cells, and tells
   the partitioner to read the load from it
   (:src:`src/amr/load_balancing/load_balancer_manager.hpp#L23-L45`). Without
   it (load balancing not active) SAMRAI counts cells.
#. **Estimators fill it.** Each level has an estimator
   (:src:`src/simulator/simulator.hpp#L275-L295`):

   * hybrid levels, ``mode="nppc"``: number of domain particles of all
     populations in each cell
     (``ConcreteLoadBalancerStrategyNPPC::compute``,
     :src:`src/amr/load_balancing/concrete_load_balancer_strategy_nppc.hpp#L77-L81`);
   * hybrid levels, ``mode="homogeneous"``, and all MHD levels: 1 per cell
     (``ConcreteLoadBalancerStrategyHomogeneous::compute``,
     :src:`src/amr/load_balancing/concrete_load_balancer_strategy_homogeneous.hpp#L63`).

   The array is refreshed at the end of every ``advanceLevel``, after a new
   level is initialized (except a regridded level 0, refreshed at the end of
   its next advance), and after a restart is loaded
   (``MultiPhysicsIntegrator``, :src:`src/amr/multiphysics_integrator.hpp#L558`,
   :src:`src/amr/multiphysics_integrator.hpp#L338-L341`,
   :src:`src/amr/multiphysics_integrator.hpp#L363-L367`). So it is up to
   date when SAMRAI regrids.
#. **Level 0 between regrids.** ``Integrator::advance`` asks
   ``_should_rebalance_now()`` and passes the answer to SAMRAI's
   ``advanceHierarchy(dt, rebalance_coarsest)``
   (:src:`src/amr/wrappers/integrator.hpp#L57-L64`). The answer is false
   unless tagging is on, level 0 is hybrid, and load balancing is active
   (:src:`src/amr/wrappers/integrator.hpp#L204-L211`). With ``every``, it is
   true every ``every`` steps (``cadence_rebalance_check``); with ``auto``,
   ``CoarsestRebalanceCheck`` decides (see *The level 0 check backs off*).

The tricks
----------

One tolerance, two jobs
^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** Two places need a "how uneven is acceptable" threshold: SAMRAI's
partitioner, and PHARE's own level 0 check. Two different definitions
would let the check ask for a rebalance that the partitioner then judges
unnecessary, forever.

**Fix.** Both use ``tol`` with SAMRAI's definition: the busiest rank may
carry up to ``tol`` above the average.
``CoarsestRebalanceCheck::is_imbalanced`` returns ``max/average - 1 > tol``
(:src:`src/amr/load_balancing/coarsest_rebalance_check.hpp#L33-L47`). The
partitioner gets the same value in every mode. With ``every=N`` there is
no level 0 check at all: ``cadence_rebalance_check`` asks for a
redistribution every ``N`` steps and only the partitioner uses ``tol``
(:src:`src/amr/wrappers/integrator.hpp#L106-L111`).

**History.** `#1308 <https://github.com/PHAREHUB/PHARE/pull/1308>`_ (load balancing cleanup) replaced an earlier
``1 - min/max`` test.

The level 0 check backs off
^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** Measuring level 0 needs the load of every rank, a collective
MPI gather. Doing it every step costs time for nothing when the load is
stable.

**Fix.** ``CoarsestRebalanceCheck`` waits ``interval`` steps between checks.
After a balanced check the interval is multiplied by
``next_rebalance_backoff_multiplier`` up to ``max_next_rebalance``; after a
rebalance it goes back to ``next_rebalance``
(:src:`src/amr/load_balancing/coarsest_rebalance_check.hpp#L49-L70`). Every
rank holds the same counters, so all ranks agree on when the gather happens
(comment at :src:`src/amr/load_balancing/coarsest_rebalance_check.hpp#L16-L17`):
none is left waiting in a collective the others skipped.

Summing the level 0 load ourselves
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** SAMRAI computes a level's load in
``CascadePartitioner::computeNonUniformWorkLoad``, which is private.

**Fix.** ``Integrator::computeNonUniformWorkLoadForLevel0`` does the same sum
with the public ``BalanceUtilities::computeNonUniformWorkload`` on each
level 0 patch, and ``mpi::collect`` gathers one number per rank
(:src:`src/amr/wrappers/integrator.hpp#L90-L119`).

Not with refinement boxes
^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** Re-partitioning with user-given refinement boxes made SAMRAI
abort in ``RefineSchedule`` ("source level finer than the destination
level").

**Fix.** ``Simulator::finalize_init`` refuses the combination with
"LoadBalancer is not supported with refinement boxes"
(:src:`src/simulator/simulator.hpp#L264-L273`), and ``pharein`` keeps load
balancing off with boxes unless the user declares one.

After a level 0 rebalance, every level is stale
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** After re-partitioning level 0, SAMRAI calls
``resetHierarchyConfiguration`` with levels (0, 0). But finer levels read
level 0 (directly or through coarse-to-fine interpolation), so their
communication :term:`schedules <schedule>` are stale too.

**Fix.** ``MultiPhysicsIntegrator::resetHierarchyConfiguration`` notifies the
messenger of *every* level, not only the range SAMRAI passes
(:src:`src/amr/multiphysics_integrator.hpp#L372-L377`). Details in
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
