.. _dev-coarsening:

==========
Coarsening
==========

.. container:: breadcrumb

   You are here: :doc:`/dev/index` › :doc:`/dev/one_time_step` › Coarsening

**In one sentence:** :term:`coarsening` copies the fine solution onto the
coarse level underneath when the two meet in time, so the coarse level
"feels" the better fine solution (:ref:`dev-one-time-step-sync`).


The rules
---------

.. plot::
   :include-source: false
   :alt: Coarsening in 2D: a coarse x-edge gets the average of its two fine x-edges; a coarse node gets the value of the fine node on top of it.

   from phare_draw import *

   fig, (ax,) = strip(1, ["one coarse cell = 2 x 2 fine cells"], phone=True)
   for v in (0, 1, 2):
       ax.plot([v, v], [0, 2], color=C["light"], lw=1)
       ax.plot([0, 2], [v, v], color=C["light"], lw=1)
   ax.plot([0, 2, 2, 0, 0], [0, 0, 2, 2, 0], color=C["ink"], lw=2.5)
   # Ex on x-edges: two fine edges along the coarse edge y = 0
   for x, v in ((0.5, 0.42), (1.5, 0.58)):
       mark(ax, "primal", x, 0, color=C["E"], size=9)
       text(ax, x, -0.28, f"{v:.2f}", color=C["E"])
   mark(ax, "dual", 1.0, 0, color=C["E"], size=13)
   text(ax, 1.0, 0.3, "coarse Ex = 0.50", color=C["E"])
   # a node: coarse node (2, 2) is also a fine node
   mark(ax, "dual", 2, 2, color=C["rho"], size=13)
   mark(ax, "primal", 2, 2, color=C["rho"], size=7)
   text(ax, 1.3, 2.45, "density node:\nfine 1.13 → coarse 1.13", color=C["rho"])
   ax.plot([-0.6, 2.6], [-0.6, 2.85], alpha=0)
   finish(fig, [ax])

Squares: coarse values, dots: fine values. Look at the edge (averaged) and the
corner node (copied).

.. list-table::
   :header-rows: 1
   :widths: 30 42 28

   * - Policy
     - Rule
     - Used for
   * - ``ElectricFieldCoarsener``
     - average of the 2 fine edges along the coarse edge; components primal
       everywhere (Ez in 2D; Ey, Ez in 1D) copy the node value
     - E, refluxed E (hybrid, MHD)
   * - ``MomentsCoarsener``
     - copy of the fine node on the coarse node (primal everywhere)
     - charge density, charge flux
   * - ``MHDFluxCoarsener``
     - average of the fine faces in the coarse face
     - MHD fluxes (density, momentum, energy)

Rules: :src:`src/amr/data/field/coarsening/electric_field_coarsener.hpp#L46-L145`,
:src:`src/amr/data/field/coarsening/moments_coarsener.hpp#L42-L78`,
:src:`src/amr/data/field/coarsening/mhd_flux_coarsener.hpp#L47-L149`; loop in
``FieldCoarsenOperator`` (:src:`src/amr/data/field/coarsening/field_coarsen_operator.hpp`).
``DefaultFieldCoarsener`` exists but is not registered.


What is coarsened
-----------------

* **Hybrid** (``synchronize``,
  :src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L626-L637`;
  :src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L912-L924`):
  E, total charge density, total charge flux. Then :term:`reflux`: the
  time-averaged fine E replaces the coarse :math:`E_\text{avg}` and coarse B
  is recomputed.
* **MHD**: ``synchronize`` is empty; ``reflux`` coarsens the time-summed fine
  fluxes and E, and the coarse level redoes its update
  (:src:`src/amr/messengers/mhd_messenger.hpp#L360-L367`,
  :src:`src/amr/messengers/mhd_messenger.hpp#L603-L615`,
  :src:`src/amr/solvers/solver_mhd.hpp#L345-L357`).


The tricks
----------

B is recomputed, not coarsened
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** Averaging fine B would not keep the coarse B the curl update of
its own E.

**Fix.** Coarse Faraday is redone with the fine E where the fine level is;
since a coarse edge is the average of its fine edges, the flux through a
coarse face changes by the sum of its fine faces' changes, and div B is
unchanged (:ref:`dev-one-time-step-sync`).

Moments are copied, not averaged
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** Averaging fine moments around a coarse node is not needed.

**Fix.** Moments are primal, so every coarse node is a fine node: its value
is copied (asserted in
:src:`src/amr/data/field/coarsening/moments_coarsener.hpp#L60-L77`).

No ghost refill for the charge flux
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** After coarsening, neighbours' ghost copies of the changed cells
are stale.

**Fix.** ``postSynchronize`` refills E and charge density ghosts, not the
charge flux: the electron velocity built from it is only used on the cells
(:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L917-L918`,
:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L652-L663`).


Go deeper
---------

* :ref:`dev-one-time-step-sync`: the order of synchronization and refluxing.
* :doc:`/dev/messengers`: synchronizers and their schedules.

Related
-------

* :doc:`/dev/field_refinement`: the opposite direction.
* :doc:`/dev/moments`: how the moments are deposited.
* :doc:`/dev/mhd_solver`: the MHD update redone at reflux.
