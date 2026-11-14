.. _dev-coarsening:

==========
Coarsening
==========

.. container:: breadcrumb

   You are here: :doc:`/dev/index` › :doc:`/dev/one_time_step` › Coarsening

**In one sentence:** :term:`coarsening` copies the fine solution onto the
coarse level underneath, at synchronization, so that the coarse level "feels"
the better fine solution where the two overlap.

It runs in ``standardLevelSynchronization``, finest level first, each time a
fine level has caught up in time with its coarser level
(:ref:`dev-one-time-step-sync`). It is the opposite of refinement
(:doc:`/dev/field_refinement`).


The picture
-----------

One coarse cell in 2D covers 2 × 2 fine cells. Each coarse value is made from
the fine values that sit on it.

.. plot::
   :include-source: false
   :alt: Coarsening in 2D: a coarse x-edge gets the average of its two fine x-edges; a coarse node gets the value of the fine node on top of it.

   import matplotlib.pyplot as plt

   fig, ax = plt.subplots(figsize=(5.2, 4.4))
   # fine grid (thin) and coarse cell (thick): coarse cell = [0,2] x [0,2]
   for v in (0, 1, 2):
       ax.plot([v, v], [0, 2], color="0.7", lw=1)
       ax.plot([0, 2], [v, v], color="0.7", lw=1)
   ax.plot([0, 2, 2, 0, 0], [0, 0, 2, 2, 0], color="k", lw=2.5)

   # Ex lives on x-edges: coarse x-edge y=0, x in [0,2]; fine x-edges [0,1] and [1,2]
   ax.plot([0.5, 1.5], [0, 0], "o", ms=11, color="#e07a1f", zorder=5)
   ax.plot([1.0], [0], "s", ms=12, color="#2a7ab9", zorder=6)
   ax.annotate("coarse Ex = average\nof the 2 fine Ex", xy=(1.0, 0.02), xytext=(0.15, 0.65),
               fontsize=12, arrowprops=dict(arrowstyle="->"))

   # nodes: coarse node (2,2) coincides with a fine node
   ax.plot([2], [2], "s", ms=12, color="#2a7ab9", zorder=6)
   ax.plot([2], [2], "o", ms=6, color="#e07a1f", zorder=7)
   ax.annotate("node value (Ez in 2D,\ndensity, flux): copied\nfrom the fine node",
               xy=(2, 2), xytext=(2.25, 1.1), fontsize=12, arrowprops=dict(arrowstyle="->"))

   ax.text(1, -0.35, "coarse cell = 2 x 2 fine cells", ha="center", fontsize=12)
   ax.set_xlim(-0.2, 4.2)
   ax.set_ylim(-0.55, 2.3)
   ax.set_aspect("equal")
   ax.axis("off")
   fig.tight_layout()

*Blue squares: coarse values; orange dots: fine values. A coarse edge value
is the average of the fine edges along it. A coarse node that is also a fine
node takes the fine value.*


How it works
------------

#. **What is coarsened (hybrid).** ``HybridHybridMessengerStrategy::synchronize``
   runs three coarsening schedules: the electric field E, the total ion charge
   density and the total ion charge flux
   (:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L626-L637`;
   registration in ``registerSyncComms``,
   :src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L912-L924`).
   Then ``reflux`` (:term:`refluxing <reflux>`) coarsens the time-averaged
   fine E onto the coarse
   :math:`E_\text{avg}`, and the coarse level recomputes B from it
   (:ref:`dev-one-time-step-sync`).

#. **What is coarsened (MHD).** The MHD messenger has an empty
   ``synchronize``; its ``reflux`` coarsens the time-summed fine fluxes of
   density, momentum and total energy, and the time-summed electric field,
   onto the coarse level, which then redoes its update with them
   (``MHDMessenger::reflux``, :src:`src/amr/messengers/mhd_messenger.hpp#L360-L367`;
   registrations at :src:`src/amr/messengers/mhd_messenger.hpp#L603-L615`
   and :src:`src/amr/messengers/mhd_messenger.hpp#L218`; ``SolverMHD::reflux``,
   :src:`src/amr/solvers/solver_mhd.hpp#L345-L357`).

#. **The operators.** ``FieldCoarsenOperator``
   (:src:`src/amr/data/field/coarsening/field_coarsen_operator.hpp`) loops
   over the coarse indices covered by the fine patch and calls a *policy* for
   each one:

   .. list-table::
      :header-rows: 1
      :widths: 30 40 30

      * - Policy
        - Rule
        - Used for
      * - ``ElectricFieldCoarsener``
        - each E component: average of the 2 fine edges along the coarse
          edge; a component that is primal in every direction (Ez in 2D,
          Ey and Ez in 1D) takes the fine value at the same node
        - E (hybrid and MHD), refluxed E
      * - ``MomentsCoarsener``
        - value of the fine node on the coarse node (moments are primal in
          every direction)
        - charge density, charge flux
      * - ``MHDFluxCoarsener``
        - each face flux: average of the fine faces in the coarse face
        - MHD fluxes of density, momentum, energy

   Rules in :src:`src/amr/data/field/coarsening/electric_field_coarsener.hpp#L46-L145`,
   :src:`src/amr/data/field/coarsening/moments_coarsener.hpp#L42-L78` and
   :src:`src/amr/data/field/coarsening/mhd_flux_coarsener.hpp#L47-L149`.
   A weighted policy, ``DefaultFieldCoarsener``, also exists but no messenger
   registers it at this commit.


The tricks
----------

B is not coarsened: it is recomputed
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** Averaging fine B onto coarse faces would give the coarse level
the fine flux, but it would not tell the coarse level how B evolved in the
cells around the fine region, and the coarse field would no longer be the
discrete curl update of its own E.

**Fix.** The coarse B is recomputed with Faraday's law from its value at the
start of the coarse step and an electric field that is the fine one where the
fine level exists (refluxing). Because a coarse edge gets the average of its
fine edges, the change of magnetic flux through a coarse face equals the sum
of the changes through the fine faces it contains: the coarse and fine
fluxes stay consistent, and :math:`\nabla\cdot B` is unchanged by the update.
The steps are on :ref:`dev-one-time-step-sync` (trick "Coarse B is
recomputed").

Moments are injected, not averaged
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** One could average fine moments around a coarse node.

**Fix.** ``MomentsCoarsener`` takes the fine value at the node that coincides
with the coarse node. The moments are deposited on primal nodes, and every
coarse node is also a fine node, so no interpolation is needed. The code
asserts that the quantity is primal in every direction
(:src:`src/amr/data/field/coarsening/moments_coarsener.hpp#L60-L77`).

No patch ghost refill for the charge flux
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** After coarsening, the coarse cells under the fine level have new
values, but the copies of those cells in neighbouring patches' ghost layers do
not.

**Fix.** ``postSynchronize`` refills the coarse patch ghosts of E and of the
charge density. The charge flux is not refilled, on purpose: the electron
velocity built from it is only computed and used on the patch cells
(comment at :src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L917-L918`;
``postSynchronize`` at
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
