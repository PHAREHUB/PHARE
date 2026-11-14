.. _dev-ghost-cells:

==========================================
Ghost cells: patch ghosts and level ghosts
==========================================

.. container:: breadcrumb

   You are here: :doc:`/dev/index` › :doc:`/dev/one_time_step` › Ghost cells

**In one sentence:** every patch stores a few layers of cells around its own
cells, the *ghost cells*, so that the solver can apply its stencils
(derivatives, particle shapes) near the patch edge without asking the
neighbours for each value.

Field values live on the nodes of the staggered (Yee) grid, at cell centers
or cell edges depending on the quantity (:doc:`/dev/gridlayout`). "Ghost
cells" means the layer of cells around the patch, with all the nodes they
contain: "ghost nodes" below are the nodes in that layer.

This page is about **field** ghost cells (B, E, moments). Particles have their
own ghost layers: see :doc:`/dev/level_ghost_particles` and
:doc:`/dev/grabbing_particles`.


The picture
-----------

.. plot::
   :include-source: false
   :alt: A level-1 patch, its ghost layer, the part copied from the neighbour patch and the part refined from level 0.

   import matplotlib.pyplot as plt
   from matplotlib.patches import Rectangle, Patch

   g = 2                       # ghost width (cells), interpolation order 1
   A = (8, 4, 8, 8)            # patch A: x0, y0, width, height (level-1 cells)
   B = (16, 4, 8, 8)           # patch B, its neighbour on level 1
   orange, blue = "#f4a259", "#7fb7e6"

   fig, ax = plt.subplots(figsize=(5.2, 3.9))

   # level 0 grid: one coarse cell = 2 x 2 fine cells
   for x in range(0, 33, 2):
       ax.plot([x, x], [0, 16], color="0.85", lw=0.8, zorder=0)
   for y in range(0, 17, 2):
       ax.plot([0, 32], [y, y], color="0.85", lw=0.8, zorder=0)

   # ghost box of patch A: level ghosts by default...
   ax.add_patch(Rectangle((A[0] - g, A[1] - g), A[2] + 2 * g, A[3] + 2 * g,
                          color=orange, alpha=0.85, zorder=1))
   # ...except where patch B's cells are: patch ghosts
   ax.add_patch(Rectangle((B[0], B[1]), g, B[3], color=blue, zorder=2))
   # patch A's own cells
   ax.add_patch(Rectangle(A[:2], A[2], A[3], color="white", zorder=2))

   # fine grid lines on the level-1 region and on A's ghost box
   for (x0, y0, w, h) in [(A[0] - g, A[1] - g, A[2] + 2 * g, A[3] + 2 * g),
                          (B[0], B[1], B[2], B[3])]:
       for x in range(x0, x0 + w + 1):
           ax.plot([x, x], [y0, y0 + h], color="0.7", lw=0.5, zorder=3)
       for y in range(y0, y0 + h + 1):
           ax.plot([x0, x0 + w], [y, y], color="0.7", lw=0.5, zorder=3)

   for P, name in [(A, "patch A"), (B, "patch B")]:
       ax.add_patch(Rectangle(P[:2], P[2], P[3], fill=False, lw=2.5,
                              color="k", zorder=5))
       ax.text(P[0] + P[2] / 2, P[1] + P[3] / 2, name, ha="center",
               va="center", fontsize=14, zorder=6)

   ax.legend(handles=[Patch(color=blue, label="patch ghosts: copied from patch B"),
                      Patch(color=orange, label="level ghosts: refined from level 0,\n"
                                                "interpolated in time")],
             loc="upper center", bbox_to_anchor=(0.5, -0.02), ncol=1,
             fontsize=13, frameon=False)
   ax.set_xlim(0, 32)
   ax.set_ylim(0, 16)
   ax.set_aspect("equal")
   ax.axis("off")
   ax.set_title("level 1 over level 0 (light grid)", fontsize=13)
   fig.tight_layout()

*Patch A on level 1 and its ghost layer (2 cells wide here). Where the ghost
layer lies on cells of another level-1 patch (B), values are copied from B:
these are* **patch ghosts**. *Where it lies outside level 1, there is no
level-1 data: values are refined from level 0 underneath and interpolated in
time. These are* **level ghosts**.


How it works
------------

#. **Size.** The ghost layer is 2 cells wide at interpolation order 1 and 4
   cells wide at orders 2 and 3 (the interpolation order is the order of the
   particle shape, see :doc:`/numerics/particles`) (``field_ghost_width``,
   :src:`src/core/models/options/hybrid_options.hpp#L24`). SAMRAI allocates
   every field on the patch box grown by this width
   (``FieldDataFactory``, :src:`src/amr/data/field/field_data_factory.hpp#L27`).

#. **One fill does both kinds.** The messenger fills B and E ghosts with one
   SAMRAI *schedule* (a precomputed list of transfers between patches, see
   :doc:`/dev/messengers`) per fill, of type ``RefinerType::GhostField``
   (:src:`src/amr/messengers/refiner.hpp#L86-L90`). It first copies values from
   the cells of the other patches of the same level (patch ghosts, including
   periodic copies across the domain), then refines from the next coarser
   level only where no patch of the level is (level ghosts).

#. **Same time, or interpolated in time.** A patch ghost is copied from a
   neighbour that is at the same time. A level ghost comes from the coarser
   level, which is already ahead in time: its value is interpolated between
   two coarse snapshots (see :ref:`dev-one-time-step-coarse-first`).

#. **Domain edges.** On a periodic direction, the ghost cells beyond the
   domain are copied from the other side, like patch ghosts. On a
   non-periodic direction they are set by the boundary conditions, called by
   SAMRAI through ``setPhysicalBoundaryConditions``
   (:src:`src/amr/data/field/refine/field_refine_patch_strategy.hpp#L101`).
   See :doc:`/numerics/domain_boundary_conditions`.


Who fills which ghosts, and when
--------------------------------

For a hybrid level (``HybridHybridMessengerStrategy``). Predictor 1,
predictor 2 and corrector are the three stages of the hybrid solver step
(:doc:`/dev/ppc_order`).

.. list-table::
   :header-rows: 1
   :widths: 18 41 41

   * - Quantity
     - How the ghosts are filled
     - When, in the solver step
   * - B
     - patch ghosts copied, level ghosts refined and interpolated in time
       (``fillMagneticGhosts``)
     - after each Faraday update (predictor 1, predictor 2, corrector) and
       after refluxing
   * - E
     - same as B (``fillElectricGhosts``)
     - after each average (for :math:`E_\text{avg}`, at mid-step time),
       after the corrector, and on the coarse level after synchronization
       (``postSynchronize``)
   * - J
     - **never filled**: Ampère computes J on the patch cells plus one ghost
       layer, from B whose ghosts are already filled
     - (none)
   * - Ion moments, per population
     - added up across patch borders, plus the level ghost particles'
       contribution at the level border
     - after each particle push, see :doc:`/dev/moments`
   * - Total charge density and charge flux
     - made equal on the border nodes shared by two patches
     - after each particle push, see :doc:`/dev/moments`

The functions are in
:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L351-L497`;
the comment on J is at
:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L833-L834`.
MHD levels fill the ghosts of their state (density, momentum, B, total energy)
and of E the same way (``MHDMessenger::fillStateGhosts``,
:src:`src/amr/messengers/mhd_messenger.hpp#L376-L394`).


The tricks
----------

NaN first, then copy, then refine only the NaNs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** A ghost node must take the neighbour's value if a neighbour patch
covers it, and a refined coarse value only otherwise. If refinement ran over
the whole ghost layer, it would overwrite the neighbour copies with
interpolated values, and two patches would disagree on the same point.

**Fix.** Three pieces work together:

.. mermaid::

   flowchart TD
     N["set every ghost node to NaN<br/>(not a number)"] --> C["copy from same-level neighbours<br/>patch ghosts get real values"]
     C --> R["refine from the coarser level<br/>write only nodes still NaN"]
     R --> L["what is left: level ghosts,<br/>now filled"]
     click R href "field_refinement.html" "Field refinement operators"

* Before each fill, the messenger sets all ghost nodes of the field to NaN
  (``setNaNsOnFieldGhosts``,
  :src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L969-L1002`;
  the comment there explains the whole mechanism).
* Fields are declared with ``fineBoundaryRepresentsVariable = false``
  (``FieldVariable``, :src:`src/amr/data/field/field_variable.hpp#L34-L39`),
  which makes SAMRAI do the same-level copy *before* the refinement.
* The refinement operator writes a fine node only if it is still NaN
  (``CompositeFieldRefiner``,
  :src:`src/amr/data/field/refine/composite_field_refiner.hpp#L171-L175`).

A side benefit: a ghost that nobody fills stays NaN and is easy to spot.

**History.** :pr:`1083` (fill B ghosts by refinement) replaced computing Faraday on
the ghost layer by filling B ghosts; :pr:`1386` (NaN ghosts with a mask) made the
NaN step cheap at scale.

Ghost fills never touch the patch's own cells
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** On a staggered (Yee) grid, some field nodes sit exactly on the
border between two patches. Both patches own a copy of such a node. A ghost
fill that also wrote "interior" values would replace a patch's own computed
value with its neighbour's.

**Fix.** The ghost refiners use a *fill pattern* (the rule that decides which
nodes a transfer may write) that removes the patch interior
from what may be written (``nonOverwriteInteriorTFfillPattern``, a
``TensorFieldFillPattern`` with ``overwrite_interior = false``,
:src:`src/amr/data/field/field_variable_fill_pattern.hpp#L129-L195`). The full
story, and what keeps the two copies of a border node equal, is on
:doc:`/dev/fill_patterns`.

J is computed, not communicated
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** Ohm's law needs J (and its Laplacian for hyper-resistivity) one
layer into the ghosts. Communicating J would add one more fill per step.

**Fix.** No refiner exists for J. Each patch computes it with Ampère on its
cells plus one ghost layer, from B whose ghosts were filled just before
(comment at :src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L833-L834`).
So B ghosts must always be filled before Ampère runs: the solver step does it
right after each Faraday update.


Go deeper
---------

* :doc:`/dev/fill_patterns`: which nodes a fill may write, and border nodes shared by two patches.
* :doc:`/dev/field_refinement`: how level-ghost values are computed from the coarse level.
* :doc:`/dev/messengers`: the schedules behind each fill, and when they are built.

Related
-------

* :doc:`/dev/one_time_step`: where in the time step each fill happens.
* :doc:`/dev/level_ghost_particles`: the particle counterpart of level ghosts.
* :doc:`/dev/moments`: how moment ghosts are summed across patch borders.
* :doc:`/data_analysis/patch_hierarchy`: ghost cells as seen in output data.
