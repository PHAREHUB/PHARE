.. _dev-field-refinement:

==========================
Field refinement operators
==========================

.. container:: breadcrumb

   You are here: :doc:`/dev/index` › :doc:`/dev/one_time_step` › :doc:`/dev/ghost_cells` › Field refinement

**In one sentence:** a field refinement operator computes fine-level values
from the coarse level; PHARE's operators are second order, and the one for B
keeps :math:`\nabla\cdot B = 0`.

They fill new levels and :term:`level ghost` cells (:doc:`/dev/ghost_cells`).
Particles are refined by :term:`particle splitting`
(:doc:`/dev/particle_splitting`).


The 1D rule
-----------

.. plot::
   :include-source: false
   :alt: The four fine sites of second-order refinement in 1D, for primal and dual quantities.

   from phare_draw import *

   fig, (a1, a2) = strip(2, titles=["primal direction (values on cell edges)",
                                    "dual direction (values at cell centres)"])
   for ax, coarse, fine in ((a1, [0, 2, 4, 6], [0, 1, 2, 3, 4, 5, 6]),
                            (a2, [1, 3, 5, 7], [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5])):
       grid(ax, (0, 0, 8 if ax is a2 else 6, 0.001), color=C["light"])
       for x in coarse:
           mark(ax, "dual", x, 1.0, color=C["L0"], size=11)
       for x in fine:
           mark(ax, "primal", x, 0.0, color=C["L1"], size=8)
       text(ax, -0.6, 1.0, "coarse", ha="right")
       text(ax, -0.6, 0.0, "fine", ha="right")
       ax.plot([-2.6, 8.3], [-0.95, 1.25], alpha=0)   # room for the labels
   arrow(a1, (2, 0.9), (2, 0.12), color=C["ink"], lw=1.4)
   text(a1, 1.75, 0.5, "copy", ha="right")
   arrow(a1, (2.1, 0.9), (2.95, 0.12), color=C["ink"], lw=1.4)
   arrow(a1, (3.9, 0.9), (3.05, 0.12), color=C["ink"], lw=1.4)
   text(a1, 3.0, -0.55, "average of the two")
   for x0 in (1, 3, 5):
       arrow(a2, (x0, 0.9), (3.5, 0.12), color=C["ink"], lw=1.4)
   text(a2, 3.5, -0.55, r"right child of cell I: $c_I + (c_{I+1}-c_{I-1})/8$")
   finish(fig, [a1, a2])

Grey squares: coarse values; blue dots: fine values. Look at which coarse
values each fine value is built from.

* Primal direction: copy on a coarse node, average halfway between two.
* Dual direction: :math:`c_I \mp (c_{I+1}-c_{I-1})/8` for the left and right
  child of coarse cell :math:`I`.
* 2D and 3D: the product of the 1D rules (``CompositeFill``,
  :src:`src/amr/data/field/refine/composite_field_refiner.hpp#L75-L87`; child
  row ``directionalProlongation``,
  :src:`src/core/data/grid/gridlayoutdefs.hpp#L153-L178`).
* Only fine values that are still NaN are written
  (:src:`src/amr/data/field/refine/composite_field_refiner.hpp#L171-L176`).
* E and scalars use this fill; B uses it as stage 1 of the scheme below
  (``makeRefineOperators_``,
  :src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L761-L780`).


B is refined in two stages to keep div B = 0
--------------------------------------------

One coarse cell of a 2D, divergence-free B, refined into 2 × 2 fine cells.
The numbers are computed with the same formulas as the code.

.. plot::
   :include-source: false
   :alt: Three panels: a coarse cell with its face B values; after stage 1 the fine cells have non-zero div B; after stage 2 every fine cell has div B = 0.

   import numpy as np

   # divergence-free coarse B from a potential A_z at coarse nodes (dx = dy = 1)
   n = 4
   X, Y = np.meshgrid(np.arange(n + 1), np.arange(n + 1), indexing="ij")
   A = 0.3 * X * Y + 0.2 * X**2 - 0.15 * Y**3 + 0.4 * np.sin(X) * np.cos(0.7 * Y)
   Bxc = A[:, 1:] - A[:, :-1]          # Bx(i, j+1/2)
   Byc = -(A[1:, :] - A[:-1, :])       # By(i+1/2, j)

   def child(c, J, right):             # dual direction: +-1/8 slope rule
       s = (c[J + 1] - c[J - 1]) / 8
       return c[J] + s if right else c[J] - s

   # stage 1 on coarse cell (1, 1): fine cells 2, 3; fine faces 2, 3, 4
   Bx = {(fx, fy): (child(Bxc[fx // 2], fy // 2, fy % 2) if fx % 2 == 0 else
                    0.5 * (child(Bxc[fx // 2], fy // 2, fy % 2)
                           + child(Bxc[fx // 2 + 1], fy // 2, fy % 2)))
         for fx in (2, 3, 4) for fy in (2, 3)}
   By = {(fx, fy): (child(Byc[:, fy // 2], fx // 2, fx % 2) if fy % 2 == 0 else
                    0.5 * (child(Byc[:, fy // 2], fx // 2, fx % 2)
                           + child(Byc[:, fy // 2 + 1], fx // 2, fx % 2)))
         for fx in (2, 3) for fy in (2, 3, 4)}
   h = 0.5
   cells = [(2, 2), (3, 2), (2, 3), (3, 3)]

   def div(bx, by, c):
       return (bx[(c[0] + 1, c[1])] - bx[c]) / h + (by[(c[0], c[1] + 1)] - by[c]) / h

   d1 = {c: div(Bx, By, c) for c in cells}

   # stage 2: correct the interior (odd index) faces, weights [3, 1]/8
   Bx2, By2 = dict(Bx), dict(By)
   toward = lambda i: -1 if i % 2 else 1
   for fy in (2, 3):
       t = toward(fy)
       Bx2[(3, fy)] += h * (3 * (d1[(3, fy)] - d1[(2, fy)])
                            + (d1[(3, fy + t)] - d1[(2, fy + t)])) / 8
   for fx in (2, 3):
       t = toward(fx)
       By2[(fx, 3)] += h * (3 * (d1[(fx, 3)] - d1[(fx, 2)])
                            + (d1[(fx + t, 3)] - d1[(fx + t, 2)])) / 8
   d2 = {c: div(Bx2, By2, c) for c in cells}

   from phare_draw import C, strip, text, finish
   from matplotlib.patches import Rectangle

   fig, axs = strip(3, titles=["coarse cell, B on its faces",
                               "stage 1: components refined alone",
                               "stage 2: interior faces fixed"], phone=True)

   def val(ax, x, y, v, vertical=False, hi=False):
       ax.text(x, y, f"{v:.3f}", ha="center", va="center", fontsize=11, color=C["B"],
               rotation=90 if vertical else 0,
               bbox=dict(fc=C["hi"] if hi else "white", ec="none", pad=0.6))

   def divtext(ax, x, y, v):
       v = 0.0 if abs(v) < 1e-12 else v
       text(ax, x, y, f"div B\n{v:+.3f}", color=C["ok"] if v == 0 else C["bad"])

   ax = axs[0]
   ax.add_patch(Rectangle((0, 0), 2, 2, fill=False, lw=2.5, ec=C["ink"]))
   val(ax, 0, 1, Bxc[1, 1], True); val(ax, 2, 1, Bxc[2, 1], True)
   val(ax, 1, 0, Byc[1, 1]); val(ax, 1, 2, Byc[1, 2])
   divtext(ax, 1, 1.25, 0.0)
   text(ax, 1, 0.55, "Bx: vertical faces\nBy: horizontal faces", fs=11, color=C["grey"])

   def fine(ax, bx, by, d, changed):
       ax.add_patch(Rectangle((0, 0), 2, 2, fill=False, lw=2.5, ec=C["ink"]))
       ax.plot([1, 1], [0, 2], color=C["grey"], lw=1)
       ax.plot([0, 2], [1, 1], color=C["grey"], lw=1)
       for (fx, fy), v in bx.items():
           val(ax, fx - 2, fy - 2 + 0.5, v, True, ("x", fx, fy) in changed)
       for (fx, fy), v in by.items():
           val(ax, fx - 2 + 0.5, fy - 2, v, False, ("y", fx, fy) in changed)
       for (cx, cy), v in d.items():
           divtext(ax, cx - 2 + 0.5, cy - 2 + 0.5, v)

   fine(axs[1], Bx, By, d1, set())
   fine(axs[2], Bx2, By2, d2, {("x", 3, 2), ("x", 3, 3), ("y", 2, 3), ("y", 3, 3)})
   finish(fig, axs, max_h=3.0)

Pink: B on faces (vertical numbers: Bx, horizontal: By). Red: div B ≠ 0,
green: div B = 0. Watch the four cells turn green in panel 3.

* **Stage 1** (``MagneticCompositeRefiner``,
  :src:`src/amr/data/field/refine/magnetic_composite_refiner.hpp#L41-L93`):
  each fine face gets the 1D rule of its own component. On a coarse face, the
  two fine values average to the coarse value, but fine cells get div B ≠ 0.
* **Stage 2** (``MagneticRefinePatchStrategy``, after each refinement,
  :src:`src/amr/data/field/refine/magnetic_refine_patch_strategy.hpp#L101-L123`):
  only faces inside the coarse cell change, by a weighted sum of the div B
  differences of the cell pairs they separate (2D: [3, 1]/8; 1D: 1/2; 3D:
  [7, 2, 2, 1]/24; ``correctFace``,
  :src:`src/amr/data/field/refine/magnetic_refine_patch_strategy.hpp#L245-L295`).
* Every fine cell ends with the div B of its coarse cell: here 0.
* Corrections read div B of the stage-1 field (``DivScratch``), not of faces
  already corrected.
* Method: Balsara, Samantaray and Subramanian, *Efficient WENO-based
  prolongation strategies for divergence-preserving vector fields*, Commun.
  Appl. Math. Comput., `doi:10.1007/s42967-021-00182-x
  <https://doi.org/10.1007/s42967-021-00182-x>`_.

**History.** :pr:`1278` (second-order AMR refinement), :pr:`1554` (readable
refiners, one touch-up).


The tricks
----------

Fill boxes are rounded out to whole coarse cells
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** Stage 2 needs all faces of a coarse cell; SAMRAI fill boxes can
cut one in half.

**Fix.** Stage 1 widens its region to whole coarse cells for each
component's normal direction (``gatherRegion_``,
:src:`src/amr/data/field/refine/magnetic_composite_refiner.hpp#L78-L93`).
Stage 2 rounds out, clips, and throws if a coarse cell is only partly covered
(``reconstructionRegion``,
:src:`src/amr/data/field/refine/magnetic_refine_patch_strategy.hpp#L165-L198`).

The order is a run-time choice
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** A template parameter would multiply the
:term:`build permutations <build permutation>`.

**Fix.** Operators call the stencil through ``IFieldRefineKernel``, chosen at
run time (:src:`src/amr/data/field/refine/field_refiner_kernel.hpp#L27-L61`).
Only order 2 exists; other values are refused
(:src:`src/amr/messengers/refinement_config.hpp`).


Go deeper
---------

* :doc:`/dev/ghost_cells`: when level ghosts are refined.
* :doc:`/dev/fill_patterns`: which nodes a refinement may write.

Related
-------

* :doc:`/dev/coarsening`: fine to coarse.
* :doc:`/dev/particle_splitting`: refining particles.
* :doc:`/dev/gridlayout`: primal and dual centerings.
