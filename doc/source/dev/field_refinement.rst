.. _dev-field-refinement:

==========================
Field refinement operators
==========================

.. container:: breadcrumb

   You are here: :doc:`/dev/index` › :doc:`/dev/one_time_step` › :doc:`/dev/ghost_cells` › Field refinement

**In one sentence:** a field refinement operator computes values on a fine
level from the coarse level underneath; PHARE's operators are second order
in space, and the one for B keeps :math:`\nabla\cdot B = 0`.

They run when a new level is created or regridded, and every time the
:term:`level ghost` cells of a fine level are filled
(:doc:`/dev/ghost_cells`). Particles are refined differently, by
:term:`particle splitting` (:doc:`/dev/particle_splitting`).


The picture
-----------

One direction, refinement ratio 2. A fine value is built from the coarse
values around it, with a rule that depends on where the fine point sits.

.. plot::
   :include-source: false
   :alt: The four fine sites of second-order refinement in 1D, for primal and dual quantities.

   import matplotlib.pyplot as plt

   fig, axs = plt.subplots(2, 1, figsize=(5.4, 4.6))

   # primal: coarse nodes at integers, fine nodes at half-integers too
   ax = axs[0]
   for I in range(4):
       ax.plot(I, 0, "s", ms=11, color="#2a7ab9")
   for f in range(0, 7):
       x = f / 2
       ax.plot(x, -0.5, "o", ms=8, color="#e07a1f" if f % 2 else "#f4a259")
   ax.annotate("", xy=(1, -0.42), xytext=(1, -0.05), arrowprops=dict(arrowstyle="->"))
   ax.text(0.92, -0.25, "copy", ha="right", va="center", fontsize=12)
   ax.annotate("", xy=(1.5, -0.42), xytext=(1.05, -0.05), arrowprops=dict(arrowstyle="->"))
   ax.annotate("", xy=(1.5, -0.42), xytext=(1.95, -0.05), arrowprops=dict(arrowstyle="->"))
   ax.text(1.5, -0.95, "midpoint\n(1/2, 1/2)", ha="center", fontsize=12)
   ax.set_title("primal direction (nodes on cell edges)", fontsize=13)

   # dual: coarse cell centres, two fine cells per coarse cell
   ax = axs[1]
   for I in range(4):
       ax.plot(I + 0.5, 0, "s", ms=11, color="#2a7ab9")
       ax.axvline(I, color="0.85", lw=1, zorder=0)
   for f in range(8):
       ax.plot(f / 2 + 0.25, -0.5, "o", ms=8, color="#e07a1f")
   ax.annotate("", xy=(1.75, -0.42), xytext=(0.55, -0.05), arrowprops=dict(arrowstyle="->"))
   ax.annotate("", xy=(1.75, -0.42), xytext=(1.5, -0.05), arrowprops=dict(arrowstyle="->"))
   ax.annotate("", xy=(1.75, -0.42), xytext=(2.45, -0.05), arrowprops=dict(arrowstyle="->"))
   ax.text(1.75, -0.95, r"right child: $c_I + (c_{I+1}-c_{I-1})/8$", ha="center", fontsize=12)
   ax.set_title("dual direction (cell centres)", fontsize=13)

   for ax in axs:
       ax.set_xlim(-0.3, 4.3)
       ax.set_ylim(-1.2, 0.35)
       ax.axis("off")
   axs[0].text(-0.25, 0, "coarse", va="center", ha="right", fontsize=12)
   axs[0].text(-0.25, -0.5, "fine", va="center", ha="right", fontsize=12)
   fig.tight_layout()

Blue squares: coarse values; orange dots: fine values. Along a primal
direction, a fine node on a coarse node is a copy, and a fine node between
two coarse nodes is their average. Along a dual direction, each coarse cell
:math:`I` splits into a left and a right child,
:math:`c_I \mp (c_{I+1}-c_{I-1})/8`: a linear profile through the coarse
cell, using its two neighbours for the slope. In 2D and 3D the rule is the
product of the 1D rules of each direction.


How it works
------------

#. **One generic fill.** ``CompositeFill`` gives each fine index a *site*
   per direction (on a coarse node, between two nodes, left child, right
   child) and applies the tensor product of the matching 1D rows
   (:src:`src/amr/data/field/refine/composite_field_refiner.hpp#L33-L39`,
   :src:`src/amr/data/field/refine/composite_field_refiner.hpp#L75-L87`). The
   child row is ``directionalProlongation``
   (:src:`src/core/data/grid/gridlayoutdefs.hpp#L153-L178`).

#. **Reach.** The fill reads one coarse cell on each side
   (``coarseStencilWidth``). SAMRAI uses this width to provide enough coarse
   data, including coarse ghost values, before refining.

#. **Only unset values are written.** The fill writes a fine value only if it
   is still NaN, so values already copied from a same-level neighbour stay
   (:src:`src/amr/data/field/refine/composite_field_refiner.hpp#L171-L176`;
   why this matters: :doc:`/dev/ghost_cells`).

#. **Which quantity uses which operator.** The hybrid messenger builds its
   operators once (``makeRefineOperators_``,
   :src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L761-L780`):
   E and scalar fields use the plain composite kernel; B uses the magnetic
   kernel, followed by the divergence touch-up described below.


The tricks
----------

B is refined in two stages to keep div B = 0
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** Each component of B lives on the faces of its own direction
(:term:`staggered grid`). Refining each component on its own gives fine
faces whose divergence is not zero, even if the coarse field is divergence
free. A non-zero :math:`\nabla\cdot B` near level borders is a source of
spurious forces.

**Fix.** Two stages, after the method of Balsara, Samantaray and Subramanian (reference below):

.. mermaid::

   flowchart TD
     S1["stage 1: refine each component<br/>with the composite fill,<br/>over whole coarse cells"] --> S2["stage 2: correct the fine faces<br/>inside each coarse cell"]
     S2 --> R["every fine cell has the divergence<br/>of its coarse cell: 0 if the coarse field is div-free"]

* **Stage 1** (``MagneticCompositeRefiner``,
  :src:`src/amr/data/field/refine/magnetic_composite_refiner.hpp#L41-L93`)
  fills every fine face of each component with the composite fill. The only
  difference from the plain kernel is the region: it is rounded out to whole
  coarse cells (``gatherRegion_``), because stage 2 needs all the faces of a
  coarse cell.
* **Stage 2** (``MagneticRefinePatchStrategy``, run by SAMRAI after each
  refinement in ``postprocessRefine``,
  :src:`src/amr/data/field/refine/magnetic_refine_patch_strategy.hpp#L101-L123`)
  adds a correction to the fine faces *inside* each coarse cell, so that all
  the fine cells (2, 4 or 8 of them) get the same divergence as their coarse
  cell. Faces that lie on a coarse face keep their stage-1 value. The
  correction of a face is a weighted sum of the divergence differences of the
  fine-cell pairs in the coarse cell, with weights 1/2 (1D), [3, 1]/8 (2D)
  and [7, 2, 2, 1]/24 (3D) (``correctFace``,
  :src:`src/amr/data/field/refine/magnetic_refine_patch_strategy.hpp#L245-L295`).

All corrections read the divergences computed from the stage-1 field (a
snapshot, ``DivScratch``): recomputing them while faces are being corrected
would mix corrected and uncorrected faces.

**History.** :pr:`1278` (second-order AMR refinement) introduced the kernels and
the touch-up; :pr:`1554` (readable refiners, one touch-up) rewrote them without
changing results.

SAMRAI's fill boxes cut coarse cells in half
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** The touch-up works on whole coarse cells. The boxes SAMRAI asks
to fill do not always cover whole coarse cells: with a temporary coarse level
and its one-cell ring, a coarse cell can be cut in half, and some of its faces
would never be written.

**Fix.** Both stages round the region out to whole coarse cells, then clip it
to the patch data (``reconstructionRegion``,
:src:`src/amr/data/field/refine/magnetic_refine_patch_strategy.hpp#L165-L198`;
rounding helpers in :src:`src/amr/data/field/refine/coarse_cell_round_out.hpp`).
The comment there notes the margin is exactly zero: rounding out needs one
more cell than the fill ring, and the smallest ghost width (2) is exactly
that. If the clip ever cut a coarse cell, the code throws rather than produce
an ill-posed correction.

The order is chosen at run time, not compiled in
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** PHARE compiles one version of the code per combination of
dimension, interpolation order and number of refined particles (a
:term:`build permutation`). Making the refinement order a template parameter
of the operators and messengers would multiply these versions.

**Fix.** The operators call the stencil through an interface,
``IFieldRefineKernel``, chosen at run time from the dictionary value
``simulation/AMR/refinement/order``
(:src:`src/amr/data/field/refine/field_refiner_kernel.hpp#L27-L61`). The
virtual call is paid once per box, not per point. Only order 2 exists at
this commit: any other value is refused (``RefinementConfig::FROM``,
:src:`src/amr/messengers/refinement_config.hpp`), and the Python option
``refinement_order`` accepts only 2 (:doc:`/user_guide/amr`).


Go deeper
---------

* :doc:`/dev/ghost_cells`: when level ghosts are filled with these operators.
* :doc:`/dev/fill_patterns`: which nodes a refinement may write.

Related
-------

* :doc:`/dev/coarsening`: the opposite direction, fine to coarse.
* :doc:`/dev/particle_splitting`: how particles are refined.
* :doc:`/dev/gridlayout`: primal and dual centerings.
* Balsara, Samantaray and Subramanian, *Efficient WENO-based
  prolongation strategies for divergence-preserving vector fields*, Commun.
  Appl. Math. Comput., `doi:10.1007/s42967-021-00182-x
  <https://doi.org/10.1007/s42967-021-00182-x>`_.
