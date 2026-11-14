==============
Load balancing
==============

**Goal:** keep every :term:`MPI` :term:`rank` about equally busy during a parallel run, so
no rank waits for a slower one.

.. contents:: On this page
   :local:
   :depth: 1


On a parallel run, PHARE splits each :term:`level` of the mesh into
:term:`patches <patch>` and gives each patch to one rank. The work of a
rank is the work of its patches: for a hybrid run, it grows with its
number of particles. As the simulation evolves, particles move and refined regions
appear, grow and vanish, so a split that was even at the start becomes
uneven. Load balancing moves patches between ranks to even it out again.

You usually do not need to do anything: with ``refinement="tagging"``, load
balancing is on by default. Use ``ph.LoadBalancer`` only
to change how and when it acts.

Minimal example
---------------

A 1D hybrid run with a dense clump of plasma in the middle and a sharp
current sheet that gets refined. Save it as ``lb.py`` and run it on 2 ranks
with ``mpirun -n 2 python3 lb.py``.

.. code-block:: python

   import numpy as np
   import pyphare.pharein as ph
   from pyphare.simulator.simulator import Simulator

   sim = ph.Simulation(
       cells=200,
       dl=0.5,
       time_step=0.005,
       final_time=0.5,
       refinement="tagging",       # load balancing is on by default with tagging
       max_nbr_levels=2,
       smallest_patch_size=20,
       diag_options={"format": "phareh5",
                     "options": {"dir": "lb_run", "mode": "overwrite"}},
   )

   def density(x):  # a dense clump: most particles sit in the middle
       return 0.2 + 2.0 * np.exp(-((x - 50.0) / 5.0) ** 2)
   def bx(x):    return 1.0
   def by(x):    return np.tanh((x - 50.0) / 1.0)   # sharp current sheet: gets refined
   def zero(x):  return 0.0
   def vth(x):   return 0.1

   ph.MaxwellianFluidModel(
       bx=bx, by=by, bz=zero,
       protons={"charge": 1, "density": density,
                "vbulkx": zero, "vbulky": zero, "vbulkz": zero,
                "vthx": vth, "vthy": vth, "vthz": vth,
                "nbr_part_per_cell": 100},
   )
   ph.ElectronModel(closure="isothermal", Te=0.1)
   ph.ElectromagDiagnostics(quantity="B", write_timestamps=[0.0, 0.5])

   ph.LoadBalancer(mode="nppc", tol=0.05)   # optional: these are the defaults

   Simulator(sim).run()

To see which rank owns which patch, read the patch attribute ``mpi_rank``:

.. code-block:: python

   from pyphare.pharesee.run import Run

   run = Run("lb_run")
   for t in (0.0, 0.5):
       B = run.GetB(t)
       for ilvl, lvl in B.levels(t).items():
           ranks = [int(p.attrs["mpi_rank"]) for p in lvl.patches]
           print(f"t={t}  level {ilvl}: patch owners {ranks}")

.. code-block:: text

   t=0.0  level 0: patch owners [0, 1]
   t=0.0  level 1: patch owners [0]
   t=0.5  level 0: patch owners [0, 1]
   t=0.5  level 1: patch owners [0, 1, 1]

At ``t=0.5`` the refined level (level 1) has been rebuilt and its patches
are now shared by both ranks.

What it does
------------

Two separate mechanisms move patches between ranks.

#. **At every regrid, on every level.** When a level is rebuilt by
   :term:`regridding` (see :doc:`/user_guide/amr`), its new patches are shared
   between ranks so that each rank gets about the same load. ``mode`` says
   how load is counted, ``tol`` how much imbalance is accepted. A level
   created for the first time is shared by number of cells; ``mode``
   applies from its next rebuild.
#. **Between regrids, on the coarsest level only.** The coarsest level
   (level 0) covers the whole domain and is never rebuilt by refinement,
   yet it usually holds most of the particles. So PHARE redistributes it
   from time to time: with ``auto`` (the default) only when a check finds it
   too uneven, with ``every=N`` every ``N`` steps without checking.
   ``auto``, ``every`` and ``on_init`` control *when* this happens.

.. plot::
   :include-source: False
   :caption: When level 0 is redistributed (red). With auto, the wait between checks (green, steps written above) doubles while level 0 stays balanced, up to 1000 steps.

   from phare_draw import *
   titles = ["auto (default): check, wait longer while balanced",
             "every=400: redistribute every 400 steps, no check"]
   fig, axes = strip(2, titles, height=1.9)
   s = 1 / 200.0
   auto = [(0, "ok"), (200, "ok"), (600, "ok"), (1400, "re"), (1600, "ok"), (1800, "ok"), (2200, "ok")]
   every = [(k, "re") for k in range(0, 2401, 400)]
   for ax, ev in zip(axes, [auto, every]):
       ax.plot([-0.3, 12.0], [0, 0], color=C["grey"], lw=1.2)
       for k in range(0, 2401, 400):
           ax.plot([k * s, k * s], [-0.08, 0.08], color=C["grey"], lw=1)
           text(ax, k * s, -0.38, str(k), fs=11)
       text(ax, 12.4, 0, "step", fs=11, ha="left")
       for st, kind in ev:
           mark(ax, "primal", st * s, 0, color=C["bad"] if kind == "re" else C["ok"], size=12)
       ax.plot([0, 13], [-0.9, 0.6], alpha=0)
   for (a, b) in zip(auto[:-1], auto[1:]):
       text(axes[0], (a[0] + b[0]) / 2 * s, 0.3, str(b[0] - a[0]), color=C["ink"], fs=11)
   text(axes[0], 9.0, -0.75, "redistributed: wait back to 200", color=C["bad"], fs=11)
   text(axes[0], 2.2, -0.75, "balanced: wait doubles", color=C["ok"], fs=11)
   finish(fig, axes)

How load is counted (``mode``)
   ``"nppc"`` (default) counts the particles in each cell: a rank's load is
   the number of particles it owns; the cost of a hybrid run grows with
   it. ``"homogeneous"`` counts each cell as 1: a rank's load is its
   number of cells.

When is a rank "too loaded" (``tol``)
   The busiest rank may carry up to ``tol`` above the average load. With
   ``tol=0.05``, a rank may carry at most 5% more than the average. Every
   redistribution aims at this. With ``auto``, the same value also decides
   *whether* level 0 is redistributed.

When level 0 is redistributed (``auto``, ``every``, ``on_init``)
   See the drawing above: ``auto`` checks at the first step (``on_init``), then less and less often while
   level 0 stays balanced; ``every=N`` redistributes every ``N`` steps without checking. Defaults below.

Options
-------

.. list-table::
   :header-rows: 1
   :widths: 25 55 20

   * - Option
     - Meaning
     - Default
   * - ``mode``
     - How load is counted: ``"nppc"`` (particles) or ``"homogeneous"``
       (cells). Applies to hybrid levels; MHD levels always count cells.
     - ``"nppc"``
   * - ``tol``
     - Accepted imbalance: the busiest rank may carry up to ``tol`` above
       the average. Target of every redistribution; with ``auto``, also
       the threshold that triggers one on level 0.
     - ``0.05``
   * - ``on_init``
     - Act on level 0 at the very first step: check it (``auto``) or
       redistribute it (``every``).
     - ``True``
   * - ``auto``
     - Choose the checking interval of level 0 automatically (see above).
       Cannot be combined with ``every``.
     - ``True`` if ``every`` is not given
   * - ``every``
     - Redistribute level 0 every ``every`` steps, without a check.
     - not set
   * - ``next_rebalance``
     - ``auto`` only: first and reset interval between level 0 checks, in
       steps.
     - ``200``
   * - ``next_rebalance_backoff_multiplier``
     - ``auto`` only: factor applied to the interval after a balanced check.
     - ``2``
   * - ``max_next_rebalance``
     - ``auto`` only: longest interval between level 0 checks, in steps.
     - ``1000``

When is load balancing active?
------------------------------

.. list-table::
   :header-rows: 1
   :widths: 45 55

   * - Your setup
     - What happens
   * - No ``refinement`` option (the default, one level)
     - Off: level 0 is shared by number of cells at start and never
       redistributed. A ``LoadBalancer`` declared here has no effect.
   * - ``refinement="tagging"``, no ``LoadBalancer``
     - On, with the defaults above.
   * - ``refinement="tagging"`` and a ``LoadBalancer``
     - On, with your values.
   * - ``refinement="boxes"`` with ``refinement_boxes``, no ``LoadBalancer``
     - Off: patches are shared between ranks by number of cells, not
       particles.
   * - ``refinement="boxes"`` with ``refinement_boxes`` and a
       ``LoadBalancer``
     - Error at start: "LoadBalancer is not supported with refinement
       boxes".
   * - MHD run (``model_options=["MHDModel"]``)
     - Mechanism 2 never runs; level 0 is redistributed only when rebuilt.

Mechanism 2 also needs ``refinement="tagging"``.

Common mistakes
---------------

* **Giving both** ``auto=True`` **and** ``every``. You get
  ``RuntimeError: LoadBalancer cannot work with both 'every' and 'auto'``.
  Pick one.
* **Declaring a** ``LoadBalancer`` **with refinement boxes.** The run stops
  at start. Remove the ``LoadBalancer`` or switch to
  ``refinement="tagging"``.
* **Expecting** ``every`` **to rebalance refined levels.** ``every``,
  ``auto`` and ``on_init`` only act on level 0. Refined levels are
  redistributed when they are rebuilt, on the refinement cadence.
* **Declaring two** ``LoadBalancer``. Only one per simulation:
  ``RuntimeError: LoadBalancer is already registered to simulation``.
* **Declaring it before** ``Simulation``. The ``LoadBalancer`` attaches to
  the current simulation: create the ``Simulation`` first, or you get
  ``RuntimeError: LoadBalancer cannot be registered as no simulation exists``.

.. seealso:: How this works inside: :doc:`/dev/load_balancing`.
