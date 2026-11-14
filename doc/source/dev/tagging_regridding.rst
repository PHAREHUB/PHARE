.. _dev-tagging-regridding:

======================
Tagging and regridding
======================

.. container:: breadcrumb

   You are here: :doc:`/dev/index` › :doc:`/dev/one_time_step` › Tagging and regridding

**In one sentence:** every second step of a level, SAMRAI rebuilds the finer
levels (:term:`regridding`); with automatic refinement, PHARE first marks the
cells where the magnetic field varies fast (:term:`tagging`), and in every
case PHARE fills the new levels with data.

PHARE writes the **tagger** and the **fill of the new level**; :term:`SAMRAI`
does the rest.

The picture
-----------

.. flow::

   SAMRAI regrid point | every 2nd step of a level
   regridAllFinerLevels
   Tag cells | PHARE applyGradientDetector, Tagger::tag
   buffer, cluster
   SAMRAI boxes | buffer, cluster, load balance | /dev/load_balancing
   initializeLevelData
   Fill the new level | messenger.regrid | /dev/messengers

Who triggers a regrid
---------------------

* ``TimeRefinementIntegrator`` sets the regrid interval from the
  :term:`refinement ratio`: 2 if even (or 1), 3 if a multiple of 3.
* ``tagCellsForRefinement`` calls PHARE's ``applyGradientDetector``: the level's
  tagger runs on every patch
  (:src:`src/amr/multiphysics_integrator.hpp#L382-L396`,
  :src:`src/amr/multiphysics_integrator.hpp#L700-L711`).
* Fixed boxes instead: ``refinement/boxes`` selects ``REFINE_BOXES``,
  ``refinement/tagging`` selects ``GRADIENT_DETECTOR``
  (:src:`src/amr/wrappers/integrator.hpp#L213-L288`).

Buffer and cluster
------------------

.. plot::
   :include-source: false
   :alt: One-dimensional schematic: two groups of tagged cells, grown by one cell, then covered by two boxes.

   from phare_draw import *
   from matplotlib.patches import Rectangle

   n = 40
   tagged = [12, 13, 30]
   buffered = sorted({c + d for c in tagged for d in (-1, 0, 1)})
   boxes_ = [(11, 14), (29, 31)]          # schematic: min and max cell of each box

   fig, axs = strip(3, ["tagged by the tagger", "grown by tag_buffer = 1", "clustered into boxes (schematic)"],
                    height=0.95, aspect=False)
   for a in axs:
       grid(a, (0, 0, n, 1))
       a.set_xlim(8, 34)
       a.set_ylim(-0.1, 1.1)
   for c in tagged:
       axs[0].add_patch(Rectangle((c, 0), 1, 1, color=C["ghost_edge"], alpha=0.9, zorder=3))
   for c in buffered:
       axs[1].add_patch(Rectangle((c, 0), 1, 1, color=C["ghost_edge"], alpha=0.35, zorder=2))
   for c in tagged:
       axs[1].add_patch(Rectangle((c, 0), 1, 1, color=C["ghost_edge"], alpha=0.9, zorder=3))
   for lo, hi in boxes_:
       axs[2].add_patch(Rectangle((lo, 0.08), hi - lo + 1, 0.84, fill=False, lw=2.5, ec=C["ink"], zorder=4))
   finish(fig, axs)

Schematic, not the output of a run: the real boxes depend on the clustering algorithm.

* Tags grow by ``tag_buffer`` cells in every direction (default 1,
  :src:`pyphare/pyphare/pharein/simulation.py#L773`).
* ``"berger"``: Berger-Rigoutsos, boxes at most 10 cells wide; ``"tile"``: SAMRAI's
  tile clustering, the default
  (:src:`src/amr/wrappers/integrator.hpp#L163-L198`,
  :src:`pyphare/pyphare/pharein/simulation.py#L597`).
* Boxes are split between ranks by the ``CascadePartitioner``
  (:src:`src/amr/wrappers/integrator.hpp#L72-L74`, :doc:`/dev/load_balancing`).

The default tagger
------------------

* ``Tagger`` is the interface (:src:`src/amr/tagging/tagger.hpp#L11-L28`);
  ``ConcreteTagger::tag`` hands the patch to a ``TaggerStrategy``
  (:src:`src/amr/tagging/concrete_tagger.hpp#L56-L67`).
* Only ``"default"`` exists (:src:`src/amr/tagging/tagger_factory.hpp#L30-L55`). To add
  a criterion, write a new ``TaggerStrategy`` subclass and select it there.
* ``DefaultTaggerStrategy::tag`` reads only B and tags a cell when the criterion exceeds
  ``tagging_threshold`` (default 0.1, :src:`pyphare/pyphare/pharein/simulation.py#L787`,
  :src:`src/amr/tagging/default_tagger_strategy.hpp#L27-L44`).

**2D and 3D.** :math:`F = B_x, B_y, B_z`, :math:`d` a direction, :math:`i` the cell along :math:`d`:

.. math::

   c_d(F) = \frac{|F_{i+2} - F_i|}{1 + |F_{i+1} - F_i|},
   \qquad
   \text{tag} \iff \max_{F,\,d}\, c_d(F) > \text{threshold}

(:src:`src/amr/tagging/default_tagger_strategy.hpp#L116-L155` in 2D,
:src:`src/amr/tagging/default_tagger_strategy.hpp#L156-L204` in 3D).

**1D.** Only :math:`B_y`, :math:`B_z`, each first averaged over 5 values:

.. math::

   \bar F_i = \tfrac15 \sum_{k=-2}^{2} F_{i+k},
   \qquad
   c(F) = \frac{|\bar F_{i+1} - \bar F_i|}{1 + |\bar F_i|},
   \qquad
   \text{tag} \iff \sqrt{c(B_y)^2 + c(B_z)^2} > \text{threshold}

(:src:`src/amr/tagging/default_tagger_strategy.hpp#L84-L108`).

.. note::

   The 1D stencil reaches 3 cells to the right: at interpolation order 1 (2 ghost
   cells, :doc:`/dev/ghost_cells`) the last cell of a patch is not computed
   (``doLastCell``, :src:`src/amr/tagging/default_tagger_strategy.hpp#L89-L91`).

What happens to the new level
-----------------------------

.. plot::
   :include-source: false
   :alt: Schematic: a new fine level takes its data from the old fine level where they overlap, and is refined from the coarser level elsewhere.

   from phare_draw import *
   from matplotlib.patches import Rectangle

   fig, (a, b, c) = strip(3, ["coarser level", "old level", "new level: source of each cell"],
                          height=1.2, aspect=False)
   for ax in (a, b, c):
       ax.set_xlim(-0.5, 20.5)
       ax.set_ylim(-0.1, 1.1)
   for x in range(0, 21, 2):              # coarse cells: twice wider (ratio 2)
       a.plot([x, x], [0, 1], color=C["grey"], lw=0.8)
   a.plot([0, 20], [0, 0], color=C["grey"], lw=0.8)
   a.plot([0, 20], [1, 1], color=C["grey"], lw=0.8)
   b.add_patch(Rectangle((2, 0), 8, 1, fill=False, lw=2.5, ec=C["L1"], zorder=4))
   grid(b, (2, 0, 10, 1), color=C["L1"])
   # new level covers 6..16: overlap with old (6..10) copied, the rest refined
   c.add_patch(Rectangle((6, 0), 4, 1, fc=C["L1"], alpha=0.35, ec="none", zorder=2))
   c.add_patch(Rectangle((10, 0), 6, 1, fc=C["ghost_edge"], alpha=0.35, ec="none", zorder=2))
   grid(c, (6, 0, 16, 1), color=C["light"])
   c.add_patch(Rectangle((6, 0), 10, 1, fill=False, lw=2.5, ec=C["L2"], zorder=4))
   text(c, 8, 0.5, "copy from\nold level", color=C["L1"])
   text(c, 13, 0.5, "refine from\ncoarser level", color=C["ghost_edge"])
   finish(fig, [a, b, c])

Schematic: the blue part is copied from the old level, the orange part refined from the coarser level.

``initializeLevelData`` is called for every level SAMRAI creates or replaces;
``isRegridding`` tells a regrid
(:src:`src/amr/multiphysics_integrator.hpp#L293-L340`).

#. Allocate the data; ``messenger.registerLevel``.
#. ``solver.onRegrid()``: before the fill above level 0, after for level 0; the hybrid solver
   clears its cached boxes and resets the ion updater
   (:src:`src/amr/solvers/solver_ppc.hpp#L115-L119`).
#. ``levelInitializer.initialize``: ``messenger.regrid`` on a regrid, else ``initLevel``
   (:src:`src/amr/level_initializer/hybrid_level_initializer.hpp#L63-L85`).
#. Rebuild what is derived: moments, level ghost contribution, J, "old" copies
   (:src:`src/amr/level_initializer/hybrid_level_initializer.hpp#L87-L169`).
#. Estimate the load (not for level 0)
   (:src:`src/amr/multiphysics_integrator.hpp#L338-L341`).

``HybridHybridMessengerStrategy::regrid``
(:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L252-L289`):

* **B**: one schedule from new, old and coarser level; a NaN check keeps copied faces
  (:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L958-L966`,
  :src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L176-L181`).
* **E**, **domain particles**: same rule, see :doc:`/dev/field_refinement`,
  :doc:`/dev/particle_splitting`.
* **Level ghost particles (old)**: ``fillTransient`` from the split of the coarser
  domain particles (above level 0); the "new" ones are refined in the next first step
  (:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L266-L282`).

Load balancing of level 0
-------------------------

``Integrator::_should_rebalance_now``: false unless tagging is on, level 0 is hybrid
and load balancing is active
(:src:`src/amr/wrappers/integrator.hpp#L205-L211`, :doc:`/dev/load_balancing`).

The tricks
----------

The skipped last cell must still be written
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** SAMRAI pre-sets to 1 the tags under the existing finer level; a skipped
cell keeps tag 1, so a refined region never goes away.

**Fix.** The skipped last cell is set to 0
(:src:`src/amr/tagging/default_tagger_strategy.hpp#L110-L114`). :pr:`1531`.

Level ghost particles are not regridded
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** Copying them from the old level led to occasional failures of the
SAMRAI MPI module, according to the code comment (issue :issue:`604`).

**Fix.** ``regrid`` uses ``fillTransient`` only
(:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L266-L281`).

Some nodes are not copied correctly by SAMRAI
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** Per the comment, a SAMRAI bug may leave wrong nodes on the ghost box
border where it touches an old-level patch (SAMRAI pull request 293).

**Fix.** None; the code only carries the warning
(:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L284-L287`).

Schedules built once per regrid
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** SAMRAI regrids the finest level first; rebuilding neighbour schedules in
``initializeLevelData`` made finer levels rebuild theirs repeatedly.

**Fix.** ``initializeLevelData`` registers only its own level; cross-level schedules
go through ``onHierarchyChanged``
(:src:`src/amr/multiphysics_integrator.hpp#L347-L378`,
:ref:`dev-messengers-schedules`). :pr:`1419`.

Go deeper
---------

* :doc:`/dev/messengers`
* :doc:`/dev/field_refinement`
* :doc:`/dev/particle_splitting`

Related
-------

* :doc:`/dev/load_balancing`
* :doc:`/dev/ghost_cells`
* :doc:`/dev/one_time_step`
