.. _dev-moments:

===================================
Moments: deposition and sum borders
===================================

.. container:: breadcrumb

   You are here: :doc:`/dev/index` › :doc:`/dev/one_time_step` ›
   :doc:`/dev/ppc_order` › Moments

**In one sentence:** the ion *moments* (density and flux of each population,
then the total charge density, bulk velocity and charge flux) are computed on
the grid nodes by *depositing* the particles of each patch, and are completed
near patch borders by summing the contributions of neighbour patches.

Why it is not trivial
---------------------

Each particle is a small cloud (its *shape*, one to a few cells wide depending on
the interpolation order). It adds a share of its weight to every node it covers.
A node near a patch border is covered by clouds of particles held by *two or
more* patches. After each patch has deposited its own particles, every copy of
that node holds only part of the total. Patches must exchange and add up these
parts. Near the border of a refined level, the plasma outside the level is also
missing; it comes from level ghost particles (:doc:`/dev/level_ghost_particles`).

This happens twice per time step, once after each of the two particle pushes of
the hybrid solver (:doc:`/dev/ppc_order`).

The picture
-----------

.. plot::
   :include-source: False
   :caption: 1D, two patches sharing the node at x = 5. After deposition, patch
      A holds only the blue parts, patch B only the orange parts, on the nodes
      both of them store (A's ghost nodes on the right, B's on the left, and the
      shared border node). The border sum gives each patch the full bar.

   import numpy as np
   import matplotlib.pyplot as plt

   fig, ax = plt.subplots(figsize=(5.5, 3.6))
   nodes = np.arange(2, 9)
   a = np.array([1.0, 1.0, 0.9, 0.55, 0.15, 0.0, 0.0])
   b = np.array([0.0, 0.0, 0.1, 0.45, 0.85, 1.0, 1.0])
   ax.bar(nodes, a, width=0.5, color="#2a7ab9", label="deposited by patch A")
   ax.bar(nodes, b, width=0.5, bottom=a, color="#f4a259", label="deposited by patch B")
   ax.axvline(5, color="k", lw=2)
   ax.text(3.3, 1.42, "patch A", ha="center", fontsize=13)
   ax.text(6.7, 1.42, "patch B", ha="center", fontsize=13)
   ax.set_ylim(0, 1.6)
   ax.set_xticks(nodes)
   ax.set_xlabel("grid nodes", fontsize=13)
   ax.set_ylabel("density", fontsize=13)
   ax.tick_params(labelsize=12)
   ax.legend(fontsize=12, frameon=False, loc="lower left", bbox_to_anchor=(0, 1.02), ncol=1)
   for s in ("top", "right"):
       ax.spines[s].set_visible(False)
   fig.tight_layout()

The full chain, run after each of the two pushes of a time step
(``SolverPPC::moveIons_``, :src:`src/amr/solvers/solver_ppc.hpp#L590-L666`):

.. mermaid::

   flowchart TD
     R["reset density and flux<br/>of each population"] --> D["deposit the patch's particles<br/>(domain + patch ghost)"]
     D --> S["sum over patch borders<br/>density and flux, per population"]
     S --> L["add level ghost particles<br/>Old and New, time-interpolated"]
     L --> T["totals: charge density,<br/>bulk velocity, charge flux"]
     T --> M["shared ghost nodes:<br/>keep the larger value"]
     click D href "grabbing_particles.html" "Patch ghost particles"
     click L href "level_ghost_particles.html" "Level ghost particles"
     click S href "fill_patterns.html" "Fill patterns: which nodes a schedule touches"
     click M href "fill_patterns.html" "Fill patterns: which nodes a schedule touches"

How it works
------------

Deposition
^^^^^^^^^^

#. At the start of a push, the density and flux of every population are set to
   zero (``resetMoments``, :src:`src/core/numerics/moments/moments.hpp#L15-L22`,
   called in ``IonUpdater::updatePopulations``,
   :src:`src/core/numerics/ion_updater/ion_updater.hpp#L96`).
#. For each particle, the ``Interpolator`` finds the first node of its cloud
   and the weights of the cloud on the nodes it covers (2, 3 or 4 nodes per
   direction at orders 1, 2, 3; weight formulas at
   :src:`src/core/numerics/interpolator/interpolator.hpp#L56-L125`).
#. It adds, on each of those nodes, ``weight × coef × w_node`` to the density
   and ``weight × coef × v × w_node`` to the three flux components
   (:src:`src/core/numerics/interpolator/interpolator.hpp#L616-L637`). ``coef``
   is 1, except for level ghost particles, which use the time weights
   :math:`1-\alpha` and :math:`\alpha`. Moments live on primal nodes, the cell
   corners (:doc:`/dev/gridlayout`).

Which particles are deposited depends on the push: see
:doc:`/dev/ppc_order`. In the real push, a patch deposits its domain particles
and its patch ghost particles, the ones that just left it for a neighbour
(:src:`src/core/numerics/ion_updater/ion_updater.hpp#L306-L309`).

Sum over patch borders
^^^^^^^^^^^^^^^^^^^^^^

#. For each population, ``fillDensityBorders`` and ``fillFluxBorders``
   (:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L388-L434`)
   first copy the patch's density (or flux) into a scratch field,
   ``sumField_`` (or ``sumVec_``).
#. A SAMRAI schedule over the level then *adds* (``PlusEqualsOp``) the
   neighbours' original values into the scratch field
   (:src:`src/amr/messengers/refiner.hpp#L100-L121`). The nodes concerned are
   the overlap of the two patches' *ghost boxes* (each patch box grown by its
   ghost cells, :doc:`/dev/ghost_cells`)
   (``FieldGhostInterpOverlapFillPattern``,
   :src:`src/amr/data/field/field_variable_fill_pattern.hpp#L327-L348`); a
   patch does not add itself (L319-L320).
#. The scratch field is copied back into the population's density (or flux).

Totals
^^^^^^

After the level ghost contribution (:doc:`/dev/level_ghost_particles`), the
solver computes the totals on every node, ghosts included
(``IonUpdater::updateIons``,
:src:`src/core/numerics/ion_updater/ion_updater.hpp#L111-L117`; ``Ions``,
:src:`src/core/data/ions/ions.hpp#L85-L178`):

* charge density :math:`\rho = \sum_s q_s n_s`;
* bulk velocity :math:`\mathbf{V} = \sum_s m_s \mathbf{F}_s / \sum_s m_s n_s`
  (:math:`\mathbf{F}_s` the flux of population :math:`s`, :math:`m_s` its mass);
* charge flux :math:`\sum_s q_s \mathbf{F}_s`, used for the electron velocity.

Finally, ``fillIonBorders``
(:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L437-L446`)
runs one more schedule on the charge density and on the charge flux: on the
nodes where two patches' ghost boxes overlap, each patch keeps the larger of its
value and its neighbour's (``SetMaxOp``, :src:`src/amr/messengers/refiner.hpp#L125-L142`;
``SetMax``, :src:`src/core/utilities/types.hpp#L546-L550`). After it, every
copy of such a node holds the same value.

Why can two copies differ at this point? One case follows from the order
above: each patch adds the level ghost particles of its *own* ghost layer, after
the border sum. Near the border of a refined level, a ghost node of patch A can
be covered by level ghost particles that only the neighbour B holds; A's copy
then misses their share, B's does not. The code does not document why the
maximum was chosen (it came with #1146).

The tricks
----------

Sum into a scratch copy, not in place
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** If the schedule added the neighbours' values directly into the
density arrays, the order of the additions would matter: once patch A has added
B's part to its node, B then adds A's value, which already contains B's part,
and counts its own part twice.

**Fix.** The sum writes into a scratch field (``sumField_``, ``sumVec_``) and
reads the neighbours' untouched density arrays; the result is copied back afterwards
(comment at :src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L394-L397`).

Leavers are deposited by the patch they leave
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** A particle that crossed into a neighbour must be counted exactly
once, although it changes owner during the step.

**Fix.** The patch it left deposits it from its patch ghost array before the
exchange; the receiving patch only gets it into its domain after the border
sum, and does not deposit it in this step (:doc:`/dev/grabbing_particles`).

Each particle is held by one patch
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** The ``+=`` border sum is only right if no particle is deposited by
two patches. When SAMRAI makes same-level patches overlap, both would hold the
particles of the overlap, and the density there came out about twice too
large.

**Fix.** Each cell has a single owner patch, and particles are only kept where
the patch owns the cell (:doc:`/dev/grabbing_particles`). The border sum itself
did not change.

**History.** #1334 (overlap single owner).

Level border nodes from Old and New
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** Near the border of a refined level, part of each node's density
comes from plasma outside the level, which no fine patch holds.

**Fix.** Level ghost particles at the start and end of the coarse step are
deposited with time weights (``fillIonPopMomentGhosts``,
:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L457-L495`);
see :doc:`/dev/level_ghost_particles`.

Go deeper
---------

* :doc:`/dev/level_ghost_particles` — the moments near a refined level's border.
* :doc:`/dev/grabbing_particles` — patch ghost particles.

Related
-------

* :doc:`/dev/ppc_order` — which particles are deposited in which push.
* :doc:`/dev/fill_patterns` — how fill patterns choose the nodes a schedule
  touches.
* :doc:`/numerics/particles` — particle shapes and interpolation orders.
