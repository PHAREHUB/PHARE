.. _dev-ghost-cells:

==========================================
Ghost cells: patch ghosts and level ghosts
==========================================

.. container:: breadcrumb

   You are here: :doc:`/dev/index` › :doc:`/dev/one_time_step` › Ghost cells

**In one sentence:** every patch stores a few layers of cells around its own,
the :term:`ghost cells <ghost cell>`, so that stencils (derivatives, particle
shapes) work at the patch edge without asking the neighbours for each value.

This page is about **field** ghosts; for particles see
:doc:`/dev/level_ghost_particles` and :doc:`/dev/grabbing_particles`.
"Ghost nodes" are the field nodes inside the ghost cells (:doc:`/dev/gridlayout`).


Patch ghosts and level ghosts
-----------------------------

.. plot::
   :include-source: false
   :alt: A level-1 patch, its ghost layer, the part copied from the neighbour patch and the part refined from level 0.

   from phare_draw import *
   from matplotlib.patches import Rectangle

   fig, (ax,) = strip(1, ["patch A on level 1, its ghost layer (2 cells)"])
   for x in range(0, 33, 2):                       # level 0: cells twice wider
       ax.plot([x, x], [0, 16], color=C["light"], lw=0.8, zorder=0)
   for y in range(0, 17, 2):
       ax.plot([0, 32], [y, y], color=C["light"], lw=0.8, zorder=0)
   patch(ax, (8, 4, 16, 12), ghosts=2, label="patch A")
   patch(ax, (16, 4, 24, 12), label="patch B")
   ax.add_patch(Rectangle((16, 4), 2, 8, fc=C["E"], alpha=0.45, ec="none", zorder=3))
   text(ax, 26.5, 14.6, "level 0", color=C["L0"])
   text(ax, 22, 0.8, "blue: patch ghosts, copied from B", color=C["E"])
   text(ax, 4.5, 0.8, "orange: level ghosts,\nrefined from level 0", color=C["ghost_edge"])
   finish(fig, [ax])

Look at patch A's ghost layer: blue where patch B is, orange outside level 1.

* **Size:** 2 cells at interpolation order 1, 4 at orders 2 and 3
  (``field_ghost_width``, :src:`src/core/models/options/hybrid_options.hpp#L24`);
  fields are allocated on the patch box grown by it
  (:src:`src/amr/data/field/field_data_factory.hpp#L27`).
* **One fill does both:** a ``GhostField`` :term:`schedule` copies from the
  level (periodic copies included), then refines from the coarser level where
  no patch is (:src:`src/amr/messengers/refiner.hpp#L86-L90`).
* **Time:** patch ghosts are copied at the same time; level ghosts are
  interpolated in time (:ref:`dev-one-time-step-coarse-first`).
* **Domain edges:** periodic → copied from the other side; otherwise boundary
  conditions (``setPhysicalBoundaryConditions``,
  :src:`src/amr/data/field/refine/field_refine_patch_strategy.hpp#L101`;
  :doc:`/numerics/domain_boundary_conditions`).


The fill order: NaN, copy, refine
---------------------------------

.. plot::
   :include-source: false
   :alt: Ghost nodes of patch A set to NaN, then the ones covered by patch B copied, then the remaining ones refined.

   from phare_draw import *

   def row(ax, ghosts):
       ax.plot([-2.5, 10.5], [0, 0], color=C["light"], lw=1, zorder=0)
       ax.axvspan(0, 8, color=C["patch"], zorder=0)
       ax.axvspan(8, 10.5, color=C["faint"], zorder=0)
       for n in range(0, 9):
           mark(ax, "primal", n, 0, color=C["ink"], size=8)
       for n, (label, col) in ghosts.items():
           mark(ax, "primal", n, 0, color=col, size=11)
           text(ax, n, 0.55, label, color=col)
       text(ax, 4, -0.6, "patch A", color=C["patch_edge"])
       text(ax, 9.3, -0.6, "patch B", color=C["grey"])
       ax.plot([-3.0, 10.8], [-0.95, 1.0], alpha=0)

   nan, cp, rf = C["bad"], C["E"], C["ghost_edge"]
   fig, axs = strip(3, ["every ghost node set to NaN",
                        "copy from patch B: its ghosts get real values",
                        "refine from level 0: only nodes still NaN"])
   row(axs[0], {-2: ("NaN", nan), -1: ("NaN", nan), 9: ("NaN", nan), 10: ("NaN", nan)})
   row(axs[1], {-2: ("NaN", nan), -1: ("NaN", nan), 9: ("copy", cp), 10: ("copy", cp)})
   row(axs[2], {-2: ("", rf), -1: ("", rf), 9: ("copy", cp), 10: ("copy", cp)})
   text(axs[2], -1.5, 0.55, "refined", color=rf)
   finish(fig, axs)

Look at the right-hand ghost nodes: once copied in step 2, step 3 leaves
them alone.

* Before each fill, ghost nodes are set to NaN (``setNaNsOnFieldGhosts``,
  :src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L969-L1002`;
  the comment explains the mechanism).
* ``fineBoundaryRepresentsVariable = false`` makes SAMRAI copy *before*
  refining (:src:`src/amr/data/field/field_variable.hpp#L34-L39`).
* The refine operator writes only NaN nodes
  (:src:`src/amr/data/field/refine/composite_field_refiner.hpp#L171-L175`).
* So neighbours never disagree, and an unfilled ghost stays NaN, easy to spot.
* :pr:`1083` (fill B ghosts by refinement), :pr:`1386` (cheap NaN masking).


Who fills which ghosts
----------------------

Hybrid level; predictor 1, predictor 2, corrector are the solver stages
(:doc:`/dev/ppc_order`).

.. list-table::
   :header-rows: 1
   :widths: 22 39 39

   * - Quantity
     - How
     - When
   * - B
     - copy + refine, in time (``fillMagneticGhosts``)
     - after each Faraday update and after refluxing
   * - E
     - same (``fillElectricGhosts``)
     - after each average (:math:`E_\text{avg}` at mid-step), after the
       corrector, after synchronization (``postSynchronize``)
   * - J
     - **never filled**: Ampère computes it on the cells plus one ghost layer
     - (none)
   * - ion moments, per population
     - summed across patch borders, plus level ghost particles
     - after each push (:doc:`/dev/moments`)
   * - total charge density, flux
     - made equal on shared border nodes
     - after each push

Functions: :src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L351-L497`.
MHD levels fill their state this way (``fillStateGhosts``,
:src:`src/amr/messengers/mhd_messenger.hpp#L376-L387`); the MHD E ghosts are
only set at physical boundaries (``NoOverlapFillPattern``,
:src:`src/amr/messengers/mhd_messenger.hpp#L523-L536`).


The tricks
----------

Ghost fills never touch the patch's own nodes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** Nodes on a patch border belong to both patches; a fill must not
replace a patch's own value.

**Fix.** Ghost refiners use a :term:`fill pattern` with
``overwrite_interior = false``
(:src:`src/amr/data/field/field_variable_fill_pattern.hpp#L129-L195`); see
:doc:`/dev/fill_patterns`.

J is computed, not communicated
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** Ohm needs J (and its Laplacian) one layer into the ghosts.

**Fix.** No refiner for J: each patch runs Ampère on its cells plus one ghost
layer, from freshly filled B
(:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L834-L835`). So
B ghosts are filled right after each Faraday update.


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
