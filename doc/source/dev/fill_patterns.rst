.. _dev-fill-patterns:

=================================
Fill patterns and overlap owners
=================================

.. container:: breadcrumb

   You are here: :doc:`/dev/index` › :doc:`/dev/one_time_step` › :doc:`/dev/ghost_cells` › Fill patterns

**In one sentence:** a *fill pattern* tells SAMRAI which nodes of a
destination patch a data transfer is allowed to write; PHARE uses its own fill
patterns so that a ghost fill never overwrites a value the patch computed
itself, and so that moments deposited near a border can be added up.

When SAMRAI copies data from patch B to patch A, it first computes the
*overlap*: the set of A's nodes that will receive values. The fill pattern is
the object that computes this overlap (a ``SAMRAI::xfer::VariableFillPattern``).


The picture
-----------

One direction, two neighbour patches on the same level, ghost width 2
(interpolation order 1). The quantity is *primal* in that direction: its
nodes sit on cell edges, so node 8 is on the border and belongs to both
patches (primal and dual centerings: :doc:`/dev/gridlayout`). Only transfers
from patch B are drawn, so A's left ghost nodes, which B does not cover, are
never written here.

.. plot::
   :include-source: false
   :alt: Which nodes of patch A each fill pattern writes, for a primal quantity in 1D.

   import matplotlib.pyplot as plt

   g = 2
   A_nodes = range(0, 9)         # patch A: cells 0..7, primal nodes 0..8
   B_nodes = range(8, 17)        # patch B: cells 8..15, primal nodes 8..16
   rows = [
       ("ghost fill\n(no overwrite)", [9, 10]),
       ("init / regrid\n(overwrite)", [8, 9, 10]),
       ("moment sum\n(ghost overlap)", [6, 7, 8, 9, 10]),
   ]
   fig, ax = plt.subplots(figsize=(5.4, 3.6))
   ax.axvspan(-0.5, 8, color="#dddddd", zorder=0)
   ax.axvspan(8, 16.5, color="#c7e0f4", zorder=0)
   ax.text(4, 3.35, "patch A", ha="center", fontsize=14)
   ax.text(12, 3.35, "patch B", ha="center", fontsize=14)
   ax.axvline(8, color="k", lw=2)
   for k, (label, written) in enumerate(rows):
       y = 2 - k
       ax.text(-3.2, y, label, ha="right", va="center", fontsize=13)
       for n in range(-g, 9 + g):          # A's nodes, ghosts included
           w = n in written
           ax.plot(n, y, "o", ms=10 if w else 7,
                   color="#e07a1f" if w else ("0.35" if 0 <= n <= 8 else "white"),
                   mec="0.35", zorder=3)
   ax.text(8, -0.75, "shared\nborder node", ha="center", va="top", fontsize=13)
   ax.set_xlim(-9.5, 17)
   ax.set_ylim(-1.6, 3.7)
   ax.set_xticks(range(-2, 11, 2))
   ax.tick_params(labelsize=13)
   ax.set_yticks([])
   for s in ("top", "right", "left"):
       ax.spines[s].set_visible(False)
   ax.set_xlabel("node index of patch A", fontsize=13)
   fig.tight_layout()

*Nodes of patch A (grey: A's own nodes 0 to 8; white: its ghost nodes).
Orange: the nodes a transfer from B may write. A ghost fill writes only A's
ghost nodes, never the shared border node 8. An init or regrid fill also
writes A's own nodes. The moment sum adds B's values on every node that is in
both patches' ghost-grown boxes (5 nodes here).*


How it works
------------

#. PHARE registers each transfer with a fill pattern, in the messengers
   (:doc:`/dev/messengers`). The patterns are in
   :src:`src/amr/data/field/field_variable_fill_pattern.hpp`.

#. SAMRAI calls the pattern's ``calculateOverlap`` (copy from a same-level
   patch) and ``computeFillBoxesOverlap`` (values made by refining coarser
   data). Both return the destination nodes to write.

#. SAMRAI always asks for the interior to be overwritten. PHARE's patterns
   ignore that request and use their own setting, fixed at compile time
   (``TensorFieldFillPattern::calculateOverlap``,
   :src:`src/amr/data/field/field_variable_fill_pattern.hpp#L140-L153`).

The patterns used in PHARE:

.. list-table::
   :header-rows: 1
   :widths: 30 30 40

   * - Pattern
     - Writes
     - Used for
   * - ``FieldFillPattern`` / ``TensorFieldFillPattern``,
       ``overwrite_interior = false``
     - ghost nodes only
     - ghost fills of B and E (patch and level ghosts), patch ghosts of the
       total charge density
   * - same, ``overwrite_interior = true``
     - ghost nodes and the patch's own nodes
     - filling B on a new level (``BalgoInit``) and after regridding
       (``BregridAlgo``)
   * - ``FieldGhostInterpOverlapFillPattern`` /
       ``TensorFieldGhostInterpOverlapFillPattern``
     - nodes inside both patches' ghost-grown boxes
     - adding up (and then equalizing) moments across patch borders, see
       :doc:`/dev/moments`
   * - ``NoOverlapFillPattern``
     - nothing from other patches
     - MHD electric field: only its ghosts outside the domain are set, by the
       boundary conditions

Registrations: hybrid in
:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L159-L201`
and :src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L846-L908`;
MHD electric field at :src:`src/amr/messengers/mhd_messenger.hpp#L523-L536`.


The tricks
----------

The shared border node is never written by a ghost fill
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** Node 8 in the picture belongs to both patches. Each patch
computes it. If a ghost fill could write it, the result would depend on the
order in which SAMRAI processes the patches.

**Fix.** With ``overwrite_interior = false``, the pattern removes the
destination's own nodes from the overlap
(``computeFillBoxesOverlap``,
:src:`src/amr/data/field/field_variable_fill_pattern.hpp#L155-L194`). Each
patch keeps its own value of the border node. The comment at the top of the
file links the SAMRAI discussion behind this choice
(:src:`src/amr/data/field/field_variable_fill_pattern.hpp#L22-L31`).

When a level is created or regridded, the new patches have no values at all
yet, so B is filled with ``overwrite_interior = true``
(``BalgoInit`` at :src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L850-L851`,
``BregridAlgo`` at :src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L180-L181`).

The setting is a type, not a value
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** SAMRAI groups the items of one algorithm into *equivalence
classes* and computes the overlaps of a class once, from its first item. If
two items differ only by a runtime setting, they can land in the same class
and one of them gets the other's overlap. PHARE hit this with the three
components of B, which have different centerings but the same data type. For
example, in 1D, Bx is primal and By is dual: By has one node fewer per patch,
so a set of nodes computed for Bx does not fit By.

**Fix.** Two things that change the overlap are made part of the C++ type:

* ``overwrite_interior`` is a template parameter of the fill pattern, so the
  two settings are two types and never share a class
  (:src:`src/amr/data/field/field_variable_fill_pattern.hpp#L22-L33`).
* each centering gets its own patch data factory type
  (``CenteredFactory``, :src:`src/amr/data/field/centered_factory.hpp#L1-L60`),
  so quantities of different centerings never share a class either. The
  comment there explains why the split is by factory type: when SAMRAI needs a
  temporary coarse level, it replaces every fill pattern and regroups the
  items, and only a split by factory type survives that.

**History.** #1492 (merge MHD schedules) introduced the split by factory type,
which lets items of different centerings share one algorithm.

Overlapping patches: one owner per cell
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** Patches of one level can overlap. When proper nesting clips a new
box to a very small size, SAMRAI grows it to the minimum patch size (at least
the ghost width) without looking at the other new boxes. In the overlap, both
patches hold the same cells. Particles there would be held twice, and their
moments counted twice.

.. plot::
   :include-source: false
   :alt: Two overlapping patches; the overlap cells belong to the patch with the smaller id.

   import matplotlib.pyplot as plt
   from matplotlib.patches import Rectangle

   fig, ax = plt.subplots(figsize=(5.2, 3.0))
   P1 = (0, 0, 8, 6)    # id 1
   P2 = (6, 2, 6, 6)    # id 2, overlaps P1 on x 6..8, y 2..6
   ax.add_patch(Rectangle(P1[:2], P1[2], P1[3], color="#7fb7e6", alpha=0.6))
   ax.add_patch(Rectangle(P2[:2], P2[2], P2[3], color="#f4a259", alpha=0.6))
   ax.add_patch(Rectangle((6, 2), 2, 4, fill=False, hatch="///", lw=0))
   for P in (P1, P2):
       ax.add_patch(Rectangle(P[:2], P[2], P[3], fill=False, lw=2, color="k"))
   ax.text(3, 3, "patch id 1", ha="center", va="center", fontsize=14)
   ax.text(10, 6, "patch id 2", ha="center", va="center", fontsize=14)
   ax.annotate("overlap:\nowned by id 1\n(smaller id)", xy=(7, 4), xytext=(13, 0.8),
               fontsize=13, arrowprops=dict(arrowstyle="->"))
   ax.set_xlim(-0.5, 19)
   ax.set_ylim(-0.5, 8.5)
   ax.set_aspect("equal")
   ax.axis("off")
   fig.tight_layout()

**Fix.** Each cell covered by several patches is owned by the patch with the
smallest ``GlobalId`` (lowest MPI rank first). Every patch computes the cells
it does *not* own from its same-level neighbours, locally, without any
communication (``makeNonOwnedBoxes`` and ``makeNonOwnedBoxesFor``,
:src:`src/amr/resources_manager/amr_utils.hpp#L244-L280`). Periodic copies of
a patch carry its id, so all ranks agree. The hybrid solver builds these boxes
once per level (``SolverPPC::make_boxes``,
:src:`src/amr/solvers/solver_ppc.hpp#L158-L171`) and only keeps domain
particles in owned cells. Without overlap, the set of non-owned cells is
empty and nothing changes. How particles move to their owner is on
:doc:`/dev/grabbing_particles`.

**History.** #1334 (single owner for overlaps).


Go deeper
---------

* :doc:`/dev/moments`: the moment sums that use the ghost-overlap pattern.
* :doc:`/dev/messengers`: which transfer uses which pattern, and when it runs.

Related
-------

* :doc:`/dev/ghost_cells`: patch ghosts and level ghosts.
* :doc:`/dev/field_refinement`: how refined values are computed, and why
  they only fill nodes still unset.
* :doc:`/dev/grabbing_particles`: particles in overlapping patches.
