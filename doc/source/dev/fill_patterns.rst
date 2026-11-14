.. _dev-fill-patterns:

=================================
Fill patterns and overlap owners
=================================

.. container:: breadcrumb

   You are here: :doc:`/dev/index` › :doc:`/dev/one_time_step` › :doc:`/dev/ghost_cells` › Fill patterns

**In one sentence:** a :term:`fill pattern` decides which nodes of a patch a
transfer may write (the *overlap*); PHARE's patterns keep a patch's own values
safe during ghost fills, and let moments be summed across borders.


Which nodes each pattern writes
-------------------------------

.. plot::
   :include-source: false
   :alt: Which nodes of patch A each fill pattern writes, for a primal quantity in 1D.

   from phare_draw import *

   def row(ax, written):
       ax.axvspan(-0.5, 8, color=C["patch"], zorder=0)
       ax.axvspan(8, 12.5, color=C["faint"], zorder=0)
       ax.plot([8, 8], [-0.4, 0.4], color=C["ink"], lw=2)
       for n in range(-2, 11):
           own = 0 <= n <= 8
           col = C["ghost_edge"] if n in written else (C["ink"] if own else C["grey"])
           mark(ax, "primal", n, 0, color=col, size=11 if n in written else 7)
       text(ax, 4, -0.65, "patch A (nodes 0..8)", color=C["patch_edge"])
       text(ax, 10.8, -0.65, "patch B", color=C["grey"])
       text(ax, 8, 0.65, "shared node 8")
       ax.plot([-2.8, 12.8], [-1.0, 1.0], alpha=0)

   fig, axs = strip(3, ["ghost fill (overwrite_interior = false)",
                        "init / regrid fill (overwrite_interior = true): node 8 too",
                        "moment sum (ghost-box overlap)"])
   row(axs[0], [9, 10])
   row(axs[1], [8, 9, 10])
   row(axs[2], [6, 7, 8, 9, 10])
   finish(fig, axs)

1D, primal quantity, ghost width 2; orange: nodes of A that a transfer from B
may write. Look at node 8: only the init/regrid and moment-sum patterns touch it.

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Pattern (:src:`src/amr/data/field/field_variable_fill_pattern.hpp`)
     - Used for
   * - ``FieldFillPattern`` / ``TensorFieldFillPattern``,
       ``overwrite_interior = false``
     - ghost fills of B and E; patch ghosts of the charge density
   * - same, ``overwrite_interior = true``
     - B on a new level (``BalgoInit``) and after regrid (``BregridAlgo``)
   * - ``FieldGhostInterpOverlapFillPattern`` (and ``Tensor...``)
     - summing, then equalizing, moments across borders (:doc:`/dev/moments`)
   * - ``NoOverlapFillPattern``
     - MHD E: only physical-boundary ghosts are set

* SAMRAI always asks to overwrite the interior; PHARE's patterns use their own
  compile-time setting instead
  (:src:`src/amr/data/field/field_variable_fill_pattern.hpp#L140-L153`).
* Registrations: hybrid
  :src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L159-L201`,
  :src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L846-L908`;
  MHD E :src:`src/amr/messengers/mhd_messenger.hpp#L523-L536`.


Overlapping patches: one owner per cell
---------------------------------------

.. plot::
   :include-source: false
   :alt: Two overlapping patches; the overlap cells belong to the patch with the smaller id.

   from phare_draw import *
   from matplotlib.patches import Rectangle

   fig, (ax,) = strip(1, ["two patches of one level that overlap"], phone=True)
   patch(ax, (0, 0, 8, 6), label="patch id 1")
   patch(ax, (6, 2, 12, 8), label="patch id 2", color=C["L2"])
   ax.add_patch(Rectangle((6, 2), 2, 4, fill=False, hatch="///", ec=C["ink"], lw=0, zorder=5))
   text(ax, 6, -1.2, "hatched: owned by id 1 (smaller id)")
   ax.plot([-0.5, 12.5], [-1.8, 9.0], alpha=0)
   finish(fig, [ax])

Look at the hatched cells: both patches hold them; only id 1 keeps particles there.

* Overlaps come from SAMRAI growing tiny boxes along periodic directions
  (PHARE turns the growth off elsewhere,
  :src:`src/amr/wrappers/hierarchy.hpp#L356-L361`).
* Owner = smallest ``GlobalId``; each patch finds its non-owned cells locally,
  no communication; periodic copies carry their patch's id
  (``makeNonOwnedBoxes``, :src:`src/amr/resources_manager/amr_utils.hpp#L244-L280`).
* The hybrid solver builds these boxes once per level
  (``SolverPPC::make_boxes``, :src:`src/amr/solvers/solver_ppc.hpp#L158-L171`)
  and keeps domain particles only in owned cells (:doc:`/dev/grabbing_particles`).
* No overlap: the set is empty, nothing changes. :pr:`1334`


The tricks
----------

The shared border node is never written by a ghost fill
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** Both patches compute node 8; if a fill wrote it, the result
would depend on SAMRAI's patch order.

**Fix.** ``overwrite_interior = false`` removes the destination's own nodes:
in ``calculateOverlap`` via the field geometry for copies
(:src:`src/amr/data/field/field_variable_fill_pattern.hpp#L140-L153`,
:src:`src/amr/data/field/field_geometry.hpp#L329-L342`), in
``computeFillBoxesOverlap`` for refined values
(:src:`src/amr/data/field/field_variable_fill_pattern.hpp#L155-L194`). New
levels use ``true``, since they have no values yet
(:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L850-L851`,
:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L180-L181`).

The setting is a type, not a value
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** SAMRAI computes the overlap of a group of items once, from the
first; Bx (primal) and By (dual) in one group would share a wrong overlap.

**Fix.** ``overwrite_interior`` is a template parameter
(:src:`src/amr/data/field/field_variable_fill_pattern.hpp#L22-L33`), and each
centering has its own factory type, ``CenteredFactory``, which survives
SAMRAI's regrouping on temporary coarse levels
(:src:`src/amr/data/field/centered_factory.hpp#L1-L60`). :pr:`1492`


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
