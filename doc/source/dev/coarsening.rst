.. _dev-coarsening:

==========
Coarsening
==========

.. container:: breadcrumb

   You are here: :doc:`/dev/index` › :doc:`/dev/refluxing` › Coarsening

**Applies to:** hybrid ✓ (E, charge density, charge flux) | MHD ✓ (only the fine flux sums, inside the reflux; its state is never coarsened)

**In one sentence:** :term:`coarsening` hands some fine quantities, averaged or
copied, to the coarse level underneath when the two meet in time
(:doc:`subcycling`); magnetic field B is never coarsened: it is recomputed
(:doc:`refluxing`).

Why
---

At a meeting time the coarse level should feel the better fine solution.
Which quantities cross, and how, depends on their centering and on the model:
that is this page.


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

The same two rules on the running example of :doc:`subcycling`, in 1D (the fine
level is the fine cells 4 to 11, under the coarse cells 2 to 5; the values are
chosen for the drawing):

.. plot::
   :include-source: false
   :alt: One dimension. A coarse node takes the value of the fine node on top of it; a coarse cell value is the average of its two fine cells.

   from phare_draw import *

   fig, (a1, a2) = strip(2, titles=["primal (nodes, Ey): copy",
                                    "dual (cells, Ex): average of the 2 fine cells"], aspect=False, height=2.5)
   fineEy = [0.30, 0.34, 0.38, 0.41, 0.40, 0.37, 0.33, 0.30, 0.28]     # fine nodes 4..12
   fineEx = [0.20, 0.30, 0.50, 0.70, 0.60, 0.40, 0.20, 0.10]           # fine cells 4..11
   for ax in (a1, a2):
       text(ax, 1.75, 1.0, "coarse", ha="right"); text(ax, 1.75, 0.0, "fine", ha="right")
       ax.plot([1.0, 6.5], [-0.6, 1.5], alpha=0)
   for k, v in enumerate(fineEy):
       x = 2 + k / 2
       mark(a1, "primal", x, 0.0, color=C["E"], size=8)
       if k % 2 == 0:
           mark(a1, "dual", x, 1.0, color=C["L0"], size=11)
           arrow(a1, (x, 0.88), (x, 0.12), color=C["grey"], lw=1.2)
           text(a1, x, 1.3, f"{v:.2f}", fs=11)
   for k, v in enumerate(fineEx):
       x = 2 + k / 2 + 0.25
       mark(a2, "primal", x, 0.0, color=C["E"], size=8)
   for c in range(4):
       v = (fineEx[2 * c] + fineEx[2 * c + 1]) / 2
       x = 2.5 + c
       mark(a2, "dual", x, 1.0, color=C["L0"], size=11)
       text(a2, x, 1.3, f"{v:.2f}", fs=11)
       arrow(a2, (x - 0.22, 0.88), (x - 0.25, 0.12), color=C["grey"], lw=1.0)
       arrow(a2, (x + 0.22, 0.88), (x + 0.25, 0.12), color=C["grey"], lw=1.0)
   finish(fig, [a1, a2], align=False)

.. claim C-co-01 C-co-02

.. list-table::
   :header-rows: 1
   :widths: 24 34 20 22

   * - Policy
     - Rule
     - Used for
     - Model
   * - ``ElectricFieldCoarsener``
     - average of the 2 fine edges along the coarse edge; components primal
       everywhere (Ez in 2D; Ey, Ez in 1D) copy the node value
     - model E, and the summed :math:`E_\text{avg}` at reflux
     - hybrid (synchronize, reflux); MHD (reflux, E)
   * - ``MomentsCoarsener``
     - copy of the fine node on the coarse node (primal everywhere)
     - charge density, charge flux
     - hybrid
   * - ``MHDFluxCoarsener``
     - 1D: copy of the fine face; 2D: mean of the 2 fine faces; 3D: mean of the 4
     - the summed fluxes of density, momentum, energy
     - MHD (reflux)

.. claim C-co-01 C-co-02 C-co-03

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
  E, total charge density, total charge flux. Then :term:`reflux`
  (:doc:`refluxing`): the time-averaged fine E replaces the coarse
  :math:`E_\text{avg}` and coarse B is recomputed. After the coarsening, the
  patch ghost cells of E and of the charge density are refilled.

.. claim C-co-06 C-rx-08

* **MHD**: the state is never coarsened: ``synchronize`` is empty. ``reflux``
  coarsens the time-summed fine fluxes and E, and the coarse level redoes its update
  (:src:`src/amr/messengers/mhd_messenger.hpp#L360-L367`,
  :src:`src/amr/messengers/mhd_messenger.hpp#L603-L615`,
  :src:`src/amr/solvers/solver_mhd.hpp#L345-L357`).


The tricks
----------

B is recomputed, not coarsened
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** The coarse B must take the fine answer and stay the result of a
curl update (that is what keeps div B = 0).

**Fix.** Coarse Faraday is redone with the fine E where the fine level is
(:doc:`refluxing`). A coarse edge is the average of its 2 fine edges (or a copy, for components
primal in every direction), so the flux
through a coarse face changes by the sum of its fine faces' changes. And div B
is unchanged whatever E is, because the divergence of a discrete curl is zero.

Moments are copied, not averaged
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** Which fine value should a coarse node get?

**Fix.** Moments are primal, so every coarse node is a fine node: its value
is copied (asserted in debug builds in
:src:`src/amr/data/field/coarsening/moments_coarsener.hpp#L60-L77`).

No ghost refill for the charge flux
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** After coarsening, neighbours' ghost copies of the changed cells
are stale.

**Fix.** ``postSynchronize`` refills E and charge density ghosts, not the
charge flux: the electron velocity built from it is only computed and used on the physical
box (no ghost cells)
(:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L917-L918`,
:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L652-L663`).


This is the last page of the chain: back to :doc:`subcycling`.

Go deeper
---------

* :doc:`refluxing`: the order of synchronization and refluxing.
* :doc:`subcycling`: when the levels meet.
* :doc:`/dev/messengers`: synchronizers and their schedules.

Related
-------

* :doc:`/dev/field_refinement`: the opposite direction.
* :doc:`/dev/moments`: how the moments are deposited.
* :doc:`/dev/mhd_solver`: the MHD update redone at reflux.
