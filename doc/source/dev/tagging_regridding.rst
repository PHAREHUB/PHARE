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

PHARE writes only two parts: the **tagger** (which cells?) and the **fill of
the new level** (which data?). The rest (when to regrid, buffering,
clustering cells into boxes) is done by :term:`SAMRAI`.

The picture
-----------

.. mermaid::

   flowchart TD
     TRI["SAMRAI TimeRefinementIntegrator<br/>atRegridPoint: every 2nd step of a level"] --> GA["SAMRAI GriddingAlgorithm<br/>regridAllFinerLevels"]
     GA --> ST["SAMRAI StandardTagAndInitialize<br/>tagCellsForRefinement"]
     ST --> AGD["PHARE MultiPhysicsIntegrator<br/>applyGradientDetector"]
     AGD --> TG["Tagger::tag<br/>default criterion on B"]
     GA --> BOX["SAMRAI: buffer the tags,<br/>cluster into boxes, load balance"]
     BOX --> ILD["PHARE MultiPhysicsIntegrator<br/>initializeLevelData<br/>(once per new level)"]
     ILD --> LI["hybrid level initializer"]
     LI --> MR["messenger.regrid"]
     LI --> MOM["moments and J<br/>recomputed on the new level"]
     MR --> FR["B and E: copy from old level<br/>or refine from coarser"]
     MR --> PS["particles: copy from old level<br/>or split from coarser"]
     click MR href "messengers.html" "Messengers"
     click FR href "field_refinement.html" "How fields are refined"
     click PS href "particle_splitting.html" "How particles are split"
     click BOX href "load_balancing.html" "Load balancing"


Who triggers a regrid
---------------------

#. **Every 2nd step of a level.** In the constructor of SAMRAI's
   ``TimeRefinementIntegrator`` (file ``algs/TimeRefinementIntegrator.cpp``),
   the regrid interval of each level is set from the
   :term:`refinement ratio` to the coarser level: 2 if the ratio is even
   (or 1), 3 if it is a multiple of 3. Its function ``atRegridPoint`` is true
   when the step number of the level is a multiple of that interval. During
   the recursive advance (:doc:`/dev/one_time_step`), the integrator then calls
   ``GriddingAlgorithm::regridAllFinerLevels``.

#. **Tag, then build boxes.** For each level to regrid, the gridding algorithm
   asks ``StandardTagAndInitialize::tagCellsForRefinement`` for tags. In
   "tagging" mode that function calls the ``applyGradientDetector`` of
   PHARE's ``MultiPhysicsIntegrator``. For every patch of the level it calls
   the tagger registered for that level
   (:src:`src/amr/multiphysics_integrator.hpp#L382-L396`; the taggers are
   registered by ``addTagger_``,
   :src:`src/amr/multiphysics_integrator.hpp#L700-L711`).

#. **Buffer, cluster.** SAMRAI then grows the tagged region by ``tag_buffer``
   cells in every direction and clusters the cells into rectangular boxes
   (``bufferTagsOnLevel`` in ``mesh/GriddingAlgorithm.cpp``). PHARE chooses
   the clustering algorithm and passes the buffer width
   (:src:`src/amr/wrappers/integrator.hpp#L163-L198`): ``"berger"`` is the
   Berger-Rigoutsos algorithm with boxes at most 10 cells wide, ``"tile"`` is
   SAMRAI's tile clustering (default in Python:
   :src:`pyphare/pyphare/pharein/simulation.py#L597`). ``tag_buffer`` defaults
   to 1 (:src:`pyphare/pyphare/pharein/simulation.py#L773`). The boxes are
   split between ranks by the ``CascadePartitioner``
   (:src:`src/amr/wrappers/integrator.hpp#L72-L74`, see
   :doc:`/dev/load_balancing`).

#. **New levels.** For every new patch, SAMRAI calls PHARE's
   ``initializeLevelData`` with the old level, so PHARE knows it is a regrid
   (next section but one).


Fixed refinement boxes instead of tagging
-----------------------------------------

The user can also give the refinement boxes, so that PHARE tags nothing.
``getUserRefinementBoxesDatabase`` builds the SAMRAI input
that selects the mode: if the input has ``refinement/boxes`` the method is
``REFINE_BOXES`` and the boxes of each level ``L<n>`` are passed to SAMRAI; else
if it has ``refinement/tagging`` the method is ``GRADIENT_DETECTOR`` (which
ends in the tagger above)
(:src:`src/amr/wrappers/integrator.hpp#L213-L288`).


The default tagger
------------------

The tagger classes are small. ``Tagger`` is the interface
(:src:`src/amr/tagging/tagger.hpp#L11-L28`). ``ConcreteTagger::tag`` gets the
patch, the model and SAMRAI's tag array, and hands them to a
``TaggerStrategy`` (:src:`src/amr/tagging/concrete_tagger.hpp#L56-L67`).
``TaggerFactory::make`` builds the strategy: only the method ``"default"``
exists, for both hybrid and MHD models
(:src:`src/amr/tagging/tagger_factory.hpp#L30-L55`). To add a criterion,
write a new ``TaggerStrategy`` subclass and select it in that factory.

``DefaultTaggerStrategy::tag`` reads only the magnetic field
(``getB``, :src:`src/amr/tagging/default_tagger_strategy.hpp#L34-L44`) and
tags a cell when a criterion exceeds ``threshold``
(``tagging_threshold`` in Python, default 0.1:
:src:`pyphare/pyphare/pharein/simulation.py#L787`,
:src:`src/amr/tagging/default_tagger_strategy.hpp#L27-L30`). The criterion
depends on the dimension.

**2D and 3D.** For each of the components :math:`F = B_x, B_y, B_z` and each
direction :math:`d`, with :math:`i` the cell index along :math:`d`:

.. math::

   c_d(F) = \frac{|F_{i+2} - F_i|}{1 + |F_{i+1} - F_i|},
   \qquad
   \text{tag} \iff \max_{F,\,d}\, c_d(F) > \text{threshold}

(:src:`src/amr/tagging/default_tagger_strategy.hpp#L116-L155` in 2D,
:src:`src/amr/tagging/default_tagger_strategy.hpp#L156-L204` in 3D).

**1D.** Only :math:`B_y` and :math:`B_z` are used. Each is first averaged over
5 values, then the criterion compares two neighbouring averages:

.. math::

   \bar F_i = \tfrac15 \sum_{k=-2}^{2} F_{i+k},
   \qquad
   c(F) = \frac{|\bar F_{i+1} - \bar F_i|}{1 + |\bar F_i|},
   \qquad
   \text{tag} \iff \sqrt{c(B_y)^2 + c(B_z)^2} > \text{threshold}

(:src:`src/amr/tagging/default_tagger_strategy.hpp#L84-L108`). The source also
holds an older 1D branch disabled by ``and false``; it never runs.

.. note::

   The 1D stencil reaches 3 cells to the right of the cell. According to
   the code comment, at interpolation order 1 (ghost layer of 2 cells, see
   :doc:`/dev/ghost_cells`) it may go beyond the last ghost node for the last
   cell of the patch, so that cell is not computed (``doLastCell`` is true only if the ghost width is
   larger than 2,
   :src:`src/amr/tagging/default_tagger_strategy.hpp#L89-L91`).


What buffering and clustering do
--------------------------------

.. plot::
   :include-source: false
   :alt: One-dimensional schematic: two groups of tagged cells, grown by one cell, then covered by two boxes.

   import matplotlib.pyplot as plt
   from matplotlib.patches import Rectangle

   n = 40
   tagged = [12, 13, 30]
   buffered = sorted({c + d for c in tagged for d in (-1, 0, 1)})
   boxes = [(11, 14), (29, 31)]          # schematic: min and max cell of each box

   fig, ax = plt.subplots(3, 1, figsize=(6.0, 3.6), sharex=True)
   for a in ax:
       for x in range(n + 1):
           a.plot([x, x], [0, 1], color="0.85", lw=0.6)
       a.set_ylim(0, 1)
       a.set_yticks([])
       for s in ("top", "right", "left"):
           a.spines[s].set_visible(False)
   for c in tagged:
       ax[0].add_patch(Rectangle((c, 0), 1, 1, color="#f4a259"))
   for c in buffered:
       ax[1].add_patch(Rectangle((c, 0), 1, 1, color="#f4a259", alpha=0.45))
   for c in tagged:
       ax[1].add_patch(Rectangle((c, 0), 1, 1, color="#f4a259"))
   for lo, hi in boxes:
       ax[2].add_patch(Rectangle((lo, 0.1), hi - lo + 1, 0.8, fill=False, lw=2.5))
   ax[0].set_title("1. tagged by the tagger", fontsize=12, loc="left")
   ax[1].set_title("2. grown by tag_buffer = 1", fontsize=12, loc="left")
   ax[2].set_title("3. clustered into boxes (schematic)", fontsize=12, loc="left")
   ax[2].set_xlabel("cell index on the coarse level", fontsize=12)
   ax[2].tick_params(labelsize=12)
   ax[2].set_xlim(0, n)
   fig.tight_layout()

*Schematic, drawn by hand, not the output of a run: the real boxes depend on
the clustering algorithm.*


What happens to the new level
-----------------------------

``MultiPhysicsIntegrator::initializeLevelData`` is the one function SAMRAI
calls for every level it creates or replaces. It knows it is regridding when
the old level is given (``isRegridding``)
(:src:`src/amr/multiphysics_integrator.hpp#L293-L340`).

#. **Allocate and register.** It allocates the data of the model, solver and
   messenger on every patch of the new level, then calls
   ``messenger.registerLevel`` (same-level :term:`schedules <schedule>`, see
   :doc:`/dev/messengers`).

#. **Tell the solver.** On a regrid of a level above 0, ``solver.onRegrid()``
   is called before the level is filled; for level 0, after. For the hybrid
   solver it clears its cached boxes and resets the ion updater
   (:src:`src/amr/solvers/solver_ppc.hpp#L115-L119`).

#. **Fill the level.** ``levelInitializer.initialize`` is called. On a regrid,
   the hybrid level initializer calls ``messenger.regrid``; on a first
   creation it calls ``initLevel`` (or fills the root level from the initial
   condition)
   (:src:`src/amr/level_initializer/hybrid_level_initializer.hpp#L63-L85`).

#. **Rebuild what is derived.** Whatever the path, the initializer then
   removes particles that sit on cells owned by a neighbouring patch,
   recomputes the ion moments, adds the contribution of the level ghost
   particles on finer levels, recomputes J with Ampere's law from B, and saves
   the "old" copies needed for time interpolation
   (:src:`src/amr/level_initializer/hybrid_level_initializer.hpp#L87-L169`).

#. **Estimate the load** of the new level (not when regridding level 0)
   (:src:`src/amr/multiphysics_integrator.hpp#L338-L340`).

The data source for the regrid is SAMRAI's rule: the new level takes data
from the **old level** where the two overlap, and is **refined from the next
coarser level** elsewhere.

``HybridHybridMessengerStrategy::regrid`` does the transfers
(:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L252-L289`):

* **B** (``magneticRegriding_``): one schedule built from the new level, the
  old level and the coarser level, with ``BregridAlgo``
  (:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L958-L966`).
  The refine operator checks for NaN before writing, so faces already copied
  from the old level are not overwritten
  (:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L176-L181`).
* **E** (``electricInitRefiners_.regrid``) and **domain particles**
  (``domainParticlesRefiners_.regrid``): the same rule, with the operators
  described in :doc:`/dev/field_refinement` and
  :doc:`/dev/particle_splitting`.
* **Level ghost particles (old)**: not regridded, filled again from the
  split of the coarser domain particles, as when a new finest level is
  created (``fillTransient``, only above level 0)
  (:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L266-L281`).
* The "new" level ghost particles are refined in the next first step
  (comment at
  :src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L283`).


Load balancing of level 0
-------------------------

Separately from regrids, ``Integrator::_should_rebalance_now`` can ask SAMRAI
to re-partition level 0 at an advance. It answers false unless tagging is on,
level 0 is hybrid and load balancing is active
(:src:`src/amr/wrappers/integrator.hpp#L205-L211`). Details:
:doc:`/dev/load_balancing`.


The tricks
----------

The skipped last cell must still be written
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** In 1D at interpolation order 1, the tagger skips the last cell of
a patch. SAMRAI pre-sets to 1 the tags of the cells covered by the existing
finer level. A cell never written keeps tag 1, so a refined region that once
covered a patch edge never goes away.

**Fix.** The last cell is explicitly set to 0 when it is skipped, in
``DefaultTaggerStrategy::tag``
(:src:`src/amr/tagging/default_tagger_strategy.hpp#L110-L114`).

**History.** :pr:`1531` (1D tagging, refined patches stick).

Level ghost particles are not regridded
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** A regrid schedule copies from the old level where they overlap.
For the old level ghost particles this led to occasional failures of the
SAMRAI MPI module, according to the code comment.

**Fix.** ``regrid`` fills them with ``fillTransient`` only, from the split of
the coarser domain particles, in
``HybridHybridMessengerStrategy::regrid``
(:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L266-L281`).

**History.** Issue :issue:`604` (cited in the code comment).

Some nodes are not copied correctly by SAMRAI
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** After the regrid schedule, some nodes may be wrong, because of a
bug in SAMRAI. Per the comment, they seem to lie on the border of the ghost
box where it touches the border of an old-level patch.

**Fix.** None in this function: the code only carries the comment, as a
warning (:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L284-L288`).

**History.** SAMRAI pull request 293 (cited in the comment).

Schedules built once per regrid
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** SAMRAI regrids the finest level first. Rebuilding the schedules
of the neighbouring levels inside ``initializeLevelData`` made the finer
levels rebuild their schedules once per coarser level regridded after them.

**Fix.** ``initializeLevelData`` registers only its own level; the
cross-level schedules are handled in ``resetHierarchyConfiguration``, through
``onHierarchyChanged`` for every level (:src:`src/amr/multiphysics_integrator.hpp#L347-L378`). See
:ref:`dev-messengers-schedules`.

**History.** :pr:`1419` (schedules once per regrid).


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
