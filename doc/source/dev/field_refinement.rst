.. _dev-field-refinement:

==========================
Field refinement operators
==========================

.. container:: breadcrumb

   You are here: :doc:`/dev/index` › :doc:`/dev/subcycling` › Field refinement

**Applies to:** hybrid ✓ | MHD ✓ (same code for both; the table below says which field uses it where)

**In one sentence:** a field refinement operator computes fine-level values
from the coarse level, at one given time; PHARE's operators are second order,
and the one for B keeps :math:`\nabla\cdot B = 0`.

Why
---

A fine level has no data of its own in three places: a new level, a regridded
patch, and its :term:`level ghost` cells (:doc:`/dev/ghost_cells`). In all three
the values come from the coarse level. This page is about *how* in space; *which
coarse time* is on :doc:`time_interpolation`. Particles are refined by
:term:`particle splitting` (:doc:`/dev/particle_splitting`).


The 1D rule, on the running example
-----------------------------------

The example used on the sub-cycling pages: level 0 has 8 cells of size 1, level 1
covers coarse cells 2 to 5, so its fine cells are 4 to 11, of size 1/2. Its
left level-ghost cells are the fine cells 2 and 3 (hybrid, interpolation order 1:
two ghost cells, :src:`src/core/models/options/hybrid_options.hpp#L24`).

.. plot::
   :include-source: false
   :alt: One dimension. Left: the two fine level-ghost values of a cell-centred quantity, built from three coarse values. Right: the two fine ghost nodes of a node-centred quantity, one copied and one averaged.

   from phare_draw import *
   from matplotlib.patches import Rectangle

   By = [1.0, 1.2, 1.6, 2.0, 2.2, 2.2, 2.0, 1.6]       # coarse cells 0..7 (dual)
   Ey = [0.0, 0.1, 0.3, 0.4, 0.4, 0.3, 0.1, 0.0, 0.0]   # coarse nodes 0..8 (primal)
   left = By[1] - (By[2] - By[0]) / 8                    # left child of coarse cell 1 = fine cell 2
   right = By[1] + (By[2] - By[0]) / 8                   # right child = fine cell 3
   assert abs(left - 1.125) < 1e-12 and abs(right - 1.275) < 1e-12
   mid = (Ey[1] + Ey[2]) / 2
   assert abs(mid - 0.2) < 1e-12

   fig, (a1, a2) = strip(2, titles=["cell-centred (dual), By: three coarse values give both children",
                                    "node-centred (primal), Ey: copy, or average of two"], aspect=False, height=2.9)
   for ax in (a1, a2):
       ax.add_patch(Rectangle((1, -0.3), 1, 1.9, fc=C["ghost"], ec="none", zorder=0))   # fine ghost region: x 1..2
       text(ax, 1.5, 1.55, "level ghosts", color=C["ghost_edge"], fs=11)
       text(ax, -0.15, 1.0, "coarse", ha="right"); text(ax, -0.15, 0.0, "fine", ha="right")
       ax.plot([-1.0, 3.3], [-0.65, 1.8], alpha=0)
   # a1: coarse cell centres 0.5, 1.5, 2.5; fine ghost cell centres 1.25 and 1.75
   for x, v in ((0.5, By[0]), (1.5, By[1]), (2.5, By[2])):
       mark(a1, "dual", x, 1.0, color=C["L0"], size=11)
       text(a1, x, 1.3, f"{v}", fs=11)
   for x, v in ((1.25, left), (1.75, right)):
       mark(a1, "primal", x, 0.0, color=C["L1"], size=9)
       text(a1, x, -0.32, f"{v:.3f}", color=C["L1"], fs=11)
   for x in (0.5, 1.5, 2.5):
       arrow(a1, (x, 0.9), (1.75, 0.1), color=C["grey"], lw=1.2)
   # a2: coarse nodes 1, 2; fine ghost nodes at 1 (copy) and 1.5 (average)
   for x, v in ((1, Ey[1]), (2, Ey[2])):
       mark(a2, "dual", x, 1.0, color=C["L0"], size=11)
       text(a2, x, 1.3, f"{v}", fs=11)
   mark(a2, "primal", 1, 0.0, color=C["L1"], size=9); text(a2, 1, -0.32, f"{Ey[1]}", color=C["L1"], fs=11)
   mark(a2, "primal", 1.5, 0.0, color=C["L1"], size=9); text(a2, 1.5, -0.32, f"{mid:.1f}", color=C["L1"], fs=11)
   arrow(a2, (1, 0.9), (1, 0.12), color=C["grey"], lw=1.4)
   text(a2, 0.82, 0.5, "copy", ha="right")
   arrow(a2, (1.05, 0.9), (1.47, 0.12), color=C["grey"], lw=1.2)
   arrow(a2, (1.95, 0.9), (1.53, 0.12), color=C["grey"], lw=1.2)
   text(a2, 1.5, 0.5, "average", ha="left", fs=11)
   finish(fig, [a1, a2], align=False)

Grey squares: coarse values; blue dots: fine level-ghost values. Look at which
coarse values each fine value is built from.

* **Cell-centred (dual) direction:** coarse cell :math:`I` gives two children,
  :math:`c_I \mp (c_{I+1}-c_{I-1})/8`. Fine cell 3 is the right child of coarse
  cell 1: :math:`1.2 + (1.6-1.0)/8 = 1.275`; fine cell 2 is its left child, 1.125.
* **Node-centred (primal) direction:** a fine node on a coarse node copies it
  (fine node 2 = coarse node 1, 0.1); halfway between two it is their average
  (fine node 3, 0.2).
* In 1D, By and Bz have no face normal, so the magnetic operator gives them exactly
  this rule; only Bx gets the second stage (below).
* The fine index :math:`f` belongs to coarse index :math:`\lfloor f/2 \rfloor`,
  also when :math:`f` is negative
  (``toCoarseIndex``, :src:`src/amr/resources_manager/amr_utils.hpp#L132-L138`).

.. claim C-fr-01 C-fr-02 C-ex-01 C-ex-02 C-ex-04 C-ex-05

* **2D and 3D:** the product of the 1D rules, one per direction
  (``CompositeFill``,
  :src:`src/amr/data/field/refine/composite_field_refiner.hpp#L88-L169`; child
  row ``directionalProlongation``,
  :src:`src/core/data/grid/gridlayoutdefs.hpp#L152-L178`). Only order 2 exists, and
  only a refinement ratio of 2 (otherwise it throws,
  :src:`src/amr/data/field/refine/composite_field_refiner.hpp#L111-L114`).
* Only fine values that are still NaN are written, so a value already there (a
  copy from a neighbour patch) is never overwritten
  (:src:`src/amr/data/field/refine/composite_field_refiner.hpp#L171-L177`).

.. claim C-fr-03

* **Averages are kept:** the two children of a cell average back to the coarse
  value (the :math:`\mp 1/8` cancel), also in 2D and 3D. For the MHD variables
  :math:`\rho`, :math:`\rho \mathbf v`, :math:`E_\text{tot}` this means a refined
  coarse cell keeps its mass, momentum and energy (for cells filled from the coarse
  level only: values already there are kept)
  (:src:`src/core/data/grid/gridlayoutdefs.hpp#L134-L151`).

.. claim C-fr-12


Who is refined, with which operator
-----------------------------------

.. list-table::
   :header-rows: 1
   :widths: 9 18 36 37

   * - Model
     - Field
     - New level, regrid
     - Level ghosts (with :doc:`time_interpolation`)
   * - hybrid
     - B
     - magnetic, two stages (below)
     - magnetic, two stages
   * - hybrid
     - Bpred
     - none
     - magnetic, two stages
   * - hybrid
     - E
     - composite rule above
     - composite rule above
   * - hybrid
     - :math:`E_\text{avg}`
     - none
     - composite rule above
   * - hybrid
     - J
     - none: Ampère on each level
     - none: computed from the refined B
   * - hybrid
     - ion density and flux
     - none: deposited from particles (:doc:`/dev/level_ghost_particles`)
     - none
   * - MHD
     - :math:`\rho`, :math:`\rho \mathbf v`, :math:`E_\text{tot}`
     - composite rule (conservative variables)
     - composite rule
   * - MHD
     - B
     - magnetic, two stages
     - magnetic, two stages
   * - MHD
     - E
     - none
     - none: patch ghosts and physical boundary only
   * - MHD
     - fluxes
     - none
     - none: patch ghosts refilled after a reflux (:doc:`refluxing`)

.. claim C-fr-04 C-fr-05 C-fr-06 C-fr-07 C-fr-08 C-fr-09 C-fr-10 C-fr-11

* Hybrid operators: ``fieldRefineOp_`` and ``EfieldRefineOp_`` use the composite
  fill, ``BInitRefineOp_`` and ``BRefineOp_`` the magnetic one
  (``makeRefineOperators_``,
  :src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L761-L780`);
  MHD: :src:`src/amr/messengers/mhd_messenger.hpp#L498-L519`.

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
* Every fine cell ends with the div B of its coarse cell (when its outer faces
  come from stage 1): here 0.
* Why the drawing is 2D: in 1D, Bx is the only component with a normal direction, and a
  divergence-free 1D Bx is uniform, so stage 2 would change nothing visible.

.. claim C-fr-13 C-fr-14 C-fr-15

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


Next: :doc:`time_interpolation`, which coarse time these values come from.

Go deeper
---------

* :doc:`subcycling`: why the coarse level is ahead of the fine one.
* :doc:`/dev/ghost_cells`: when level ghosts are refined.
* :doc:`/dev/fill_patterns`: which nodes a refinement may write.

Related
-------

* :doc:`/dev/coarsening`: fine to coarse.
* :doc:`/dev/particle_splitting`: refining particles.
* :doc:`/dev/gridlayout`: primal and dual centerings.
