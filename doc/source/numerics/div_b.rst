===============================================
Divergence of B: on one level and across levels
===============================================

**Goal:** know when :math:`\nabla\cdot\mathbf B` is zero in a PHARE run, in what sense,
on every level, and why. Hybrid and MHD are both covered; each statement says which.

**In one sentence:** a step never creates or removes the divergence of B in a cell that
it advances, so the divergence of a run is set by the initial field and by the way a fine
level is created from a coarse one; refinement is built to copy the divergence of the
coarse cell, and refluxing keeps the fluxes of the two levels equal at their border.

.. contents:: On this page
   :local:
   :depth: 1

Is div B zero in my run?
------------------------

* **Initial field.** It is whatever the sampled formulas give. Nothing corrects it
  (:ref:`below <divb-initial>`).
* **One level, every step, hybrid and MHD.** The divergence of every advanced cell is
  unchanged, to round-off (:ref:`below <divb-one-level>`).
* **A new or regridded fine level.** Each fine cell gets the divergence of its coarse cell:
  zero if that is zero (:ref:`below <divb-refinement>`).
* **Ghost cells.** They are not advanced. Next to the coarse–fine border they are not
  divergence-free in general; they never change the divergence of an advanced cell.
* **Synchronization of the levels.** Each coarse cell keeps its divergence, and the magnetic
  flux through a coarse face stays equal to the sum of the fluxes through the fine faces
  inside it, for a fixed time step (:ref:`below <divb-sync>`).
* **Measured** on a two-level MHD run: :ref:`below <divb-measured>`.

Why it matters
--------------

Maxwell's equations say :math:`\nabla\cdot\mathbf B=0`. A scheme that lets it drift gives
a field with magnetic charges, and the force on the plasma then has a spurious part. The
point of this page is that PHARE's updates are built so that nothing in a time step makes
it drift.

On the Yee grid (:doc:`grid`) each component of B sits on the faces perpendicular to it,
and the discrete divergence of a cell is the sum of the differences across its faces. In
2D:

.. math::

   (\nabla\cdot\mathbf B)_{ij} = \frac{B_x^{\,i+\frac12} - B_x^{\,i-\frac12}}{\Delta x}
   + \frac{B_y^{\,j+\frac12} - B_y^{\,j-\frac12}}{\Delta y}

where :math:`B_x^{\,i\pm\frac12}` are the values of :math:`B_x` on the two faces of the cell
along x, and likewise for :math:`B_y`

(the operator is :src:`src/core/numerics/differential_operators.hpp#L159-L182`).

.. _divb-one-level:

One level: why a step keeps div B
---------------------------------

Both models advance B with Faraday's law, :math:`\mathbf B \leftarrow \mathbf B-\Delta t\,\nabla\times\mathbf E`,
with E on the cell edges (:src:`src/core/numerics/faraday/faraday.hpp#L44-L57`):

* **Hybrid:** the three Faraday updates of the predictor–predictor–corrector step all
  start from :math:`\mathbf B^n` and subtract a curl
  (:src:`src/amr/solvers/solver_ppc.hpp#L417-L423`,
  :src:`src/amr/solvers/solver_ppc.hpp#L456-L461`,
  :src:`src/amr/solvers/solver_ppc.hpp#L498-L503`).
* **MHD:** each Euler stage updates B from the edge field of that stage, and every
  integrator adds the edge fields of its stages with the same weights as the fluxes, so the
  update is the curl of a weighted E (:doc:`mhd_time_integration`;
  :src:`src/amr/solvers/time_integrator/euler_using_computed_flux.hpp#L30`,
  :src:`src/amr/solvers/time_integrator/base_mhd_timestepper.hpp#L115-L124`).

The discrete curl and the discrete divergence are made of the same one-cell differences
(:src:`src/core/numerics/differential_operators.hpp#L45-L139`), so the divergence of a
curl is zero term by term. The figure shows it on one cell.

.. plot::
   :include-source: False
   :alt: One 2D cell with E_z on its four corners; the change of B on its four faces; the change of divergence is zero.

   import matplotlib.pyplot as plt
   from matplotlib.patches import Rectangle
   from phare_draw import C

   dt, dx = 0.1, 1.0
   E = {"BL": 0.30, "BR": 0.10, "TL": -0.20, "TR": 0.50}      # E_z on the four corners
   dBx_L = -dt * (E["TL"] - E["BL"]) / dx                      # Faraday: dBx = -dt dEz/dy
   dBx_R = -dt * (E["TR"] - E["BR"]) / dx
   dBy_B = dt * (E["BR"] - E["BL"]) / dx                       # dBy = +dt dEz/dx
   dBy_T = dt * (E["TR"] - E["TL"]) / dx
   ddiv = (dBx_R - dBx_L) / dx + (dBy_T - dBy_B) / dx
   assert abs(ddiv) < 1e-12
   fig, ax = plt.subplots(figsize=(3.4, 3.4))
   ax.add_patch(Rectangle((0, 0), 1, 1, fill=False, lw=2.5, ec=C["ink"]))
   for (x, y), k in [((0, 0), "BL"), ((1, 0), "BR"), ((0, 1), "TL"), ((1, 1), "TR")]:
       ax.plot([x], [y], "o", color=C["E"], ms=8)
       ax.text(x + (-0.06 if x == 0 else 0.06), y + (-0.1 if y == 0 else 0.1), f"E_z = {E[k]:+.2f}",
               color=C["E"], fontsize=8, ha="right" if x == 0 else "left", va="center")
   for (x, y, v, rot) in [(0, 0.5, dBx_L, 90), (1, 0.5, dBx_R, 90), (0.5, 0, dBy_B, 0), (0.5, 1, dBy_T, 0)]:
       ax.text(x, y, f"ΔB {v:+.2f}", color=C["B"], fontsize=8, ha="center", va="center", rotation=rot,
               bbox=dict(fc="white", ec="none", pad=1.0))
   ax.text(0.5, 0.5, f"change of div B\n= {abs(ddiv):.3f}", ha="center", va="center", fontsize=9, color=C["ok"])
   ax.set_xlim(-0.75, 1.75); ax.set_ylim(-0.3, 1.3); ax.set_aspect("equal"); ax.set_axis_off()
   fig.tight_layout()

The four corner values are arbitrary: take any numbers, the four face changes always
cancel in the cell (:math:`\Delta t=0.1`, :math:`\Delta x=\Delta y=1` in the figure).
This is why:

* a step keeps the divergence of every cell it advances, whatever E is, and so
  whatever the ghost cells hold: ghost values only enter E;
* a step does not remove a divergence that is already there;
* in 1D nothing ever changes :math:`B_x`: the derivative along a missing direction is zero,
  so the curl has no x component (hybrid and MHD,
  :src:`src/core/numerics/differential_operators.hpp#L84-L101`), and the MHD flux of
  :math:`B_x` is zero (:src:`src/core/numerics/MHD_equations/MHD_equations.hpp#L38`).
  :math:`\nabla\cdot\mathbf B=\partial_x B_x` stays as initialised: zero only for a uniform ``bx``.

.. _divb-initial:

The initial field
-----------------

The initial B is the user function sampled at the position of each face
(:src:`src/core/data/field/initializers/field_user_initializer.hpp#L17-L40`). Nothing
corrects it afterwards, and at this commit PHARE has no vector-potential input, no
divergence-cleaning step and no run-time check of the divergence.

* If each component depends only on the other coordinates (:math:`B_x(y)`, :math:`B_y(x)`,
  as in the Orszag–Tang set-up), the discrete divergence is zero exactly.
* If the field is the **pointwise sampled curl of a potential**, it is divergence-free
  analytically but not on the grid: its discrete divergence falls like :math:`\Delta x^2`.

The next figure is a documentation run (:doc:`/validation/mhd`): a Gaussian-potential field
sampled this way on four grids, one level, 2D MHD. Dots: maximum :math:`|\nabla\cdot\mathbf B|` at
:math:`t=0`; crosses: after ten steps.

.. plot::
   :include-source: False
   :alt: Maximum divergence of B of a sampled curl field against the cell size, at t=0 and after ten steps, with a slope-2 line.

   import os
   import numpy as np
   import matplotlib.pyplot as plt
   from phare_draw import C

   for base in ("../validation/data", "validation/data", "doc/source/validation/data"):
       f = os.path.join(base, "mhd_divb_init.txt")
       if os.path.exists(f):
           break
   d = np.loadtxt(f, comments="#")
   fig, ax = plt.subplots(figsize=(3.4, 3.0))
   ax.loglog(d[:, 1], d[:, 2], "o-", color=C["L1"], label="t = 0", ms=6, lw=1.8)
   ax.loglog(d[:, 1], d[:, 3], "x", color=C["bad"], label="after 10 steps", ms=9, mew=2)
   xs = np.array([d[-1, 1], d[0, 1]])
   ax.loglog(xs, d[0, 2] * (xs / d[0, 1]) ** 2, "--", color=C["grey"], lw=1.2)
   ax.text(xs[0] * 1.1, d[0, 2] * (xs[0] / d[0, 1]) ** 2 * 1.3, "slope 2", fontsize=8, color=C["grey"])
   ax.set_xlabel("cell size Δx", fontsize=9); ax.set_ylabel("max |div B|", fontsize=9)
   ax.legend(fontsize=8, frameon=False)
   for s in ("top", "right"):
       ax.spines[s].set_visible(False)
   fig.tight_layout()

The divergence is the same after ten steps, to the precision of the output files (single
precision, which resolves a divergence to about :math:`10^{-5}` on the finest grid), and it falls by a factor of 3.5, 3.9
and 4.0 when the cell size is halved. To get less divergence initially, use a field that is
exactly divergence-free on the grid, or a finer grid.

``Run.GetDivB`` computes the divergence from the diagnostics. It works in 2D only
(:src:`pyphare/pyphare/pharesee/run/utils.py#L124-L146`).

.. admonition:: Physical boundaries (MHD)
   :class: note

   At physical boundaries the MHD boundary conditions fill the ghost B with
   divergence-free transverse conditions (the reflective, open and inflow types; they are
   refused for hybrid): see :doc:`domain_boundary_conditions`.

.. _divb-refinement:

Across levels, one coarse step
------------------------------

Take two levels with refinement ratio 2. A coarse step has these moments, each with its own
answer for div B.

.. flow::

   Create or regrid the fine level | fine B refined from the coarse level | /dev/field_refinement
   then
   Each fine sub-step | level ghosts refilled; interior advanced by Faraday | /dev/time_interpolation
   then
   End of the coarse step | the coarse level redoes its B update with the fine electric field | /dev/refluxing

Creating the fine level
~~~~~~~~~~~~~~~~~~~~~~~

A new fine level, or one rebuilt by a regrid, is *not* sampled again from the user
functions: its B is refined from the coarse level, so each fine cell starts with the
divergence of its coarse parent (hybrid:
:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L812-L831`; MHD:
:src:`src/amr/messengers/mhd_messenger.hpp#L493-L514`). Only the second-order operator exists
(``refinement_order`` must be 2).

The B refinement has two stages (hybrid and MHD):

#. **Each component alone**, with the rules of :doc:`amr`. On a coarse face, the two fine
   faces average to the coarse value. The fine cells then have a divergence that differs
   from their parent's
   (:src:`src/amr/data/field/refine/magnetic_composite_refiner.hpp#L41-L93`).
#. **A correction of the faces inside the coarse cell**, a weighted sum of the
   divergence differences of the neighbouring fine cells (weights 1/2 in 1D, 3/1 over 8 in
   2D, 7/2/2/1 over 24 in 3D), so that all fine cells of a coarse cell get the same
   divergence, equal to the coarse cell's
   (:src:`src/amr/data/field/refine/magnetic_refine_patch_strategy.hpp#L100-L124`,
   :src:`src/amr/data/field/refine/magnetic_refine_patch_strategy.hpp#L256-L300`). The faces
   on the coarse faces keep their stage-1 value.

.. plot::
   :include-source: False
   :alt: A coarse cell whose divergence is 0.3; after stage 1 its four fine cells have 0.349, 0.254, 0.346, 0.251; after stage 2 all have 0.300.

   import numpy as np
   # coarse B from a potential plus a term that gives every coarse cell the divergence q0 = 0.3
   n = 4
   X, Y = np.meshgrid(np.arange(n + 1), np.arange(n + 1), indexing="ij")
   A = 0.3 * X * Y + 0.2 * X**2 - 0.15 * Y**3 + 0.4 * np.sin(X) * np.cos(0.7 * Y)
   q0 = 0.3
   Bxc = (A[:, 1:] - A[:, :-1]) + q0 * np.arange(n + 1)[:, None] * np.ones((1, n))
   Byc = -(A[1:, :] - A[:-1, :])
   def child(c, J, right):
       s = (c[J + 1] - c[J - 1]) / 8
       return c[J] + s if right else c[J] - s
   Bx = {(fx, fy): (child(Bxc[fx // 2], fy // 2, fy % 2) if fx % 2 == 0 else
                    0.5 * (child(Bxc[fx // 2], fy // 2, fy % 2) + child(Bxc[fx // 2 + 1], fy // 2, fy % 2)))
         for fx in (2, 3, 4) for fy in (2, 3)}
   By = {(fx, fy): (child(Byc[:, fy // 2], fx // 2, fx % 2) if fy % 2 == 0 else
                    0.5 * (child(Byc[:, fy // 2], fx // 2, fx % 2) + child(Byc[:, fy // 2 + 1], fx // 2, fx % 2)))
         for fx in (2, 3) for fy in (2, 3, 4)}
   h = 0.5
   cells = [(2, 2), (3, 2), (2, 3), (3, 3)]
   def div(bx, by, c):
       return (bx[(c[0] + 1, c[1])] - bx[c]) / h + (by[(c[0], c[1] + 1)] - by[c]) / h
   d1 = {c: div(Bx, By, c) for c in cells}
   Bx2, By2 = dict(Bx), dict(By)
   toward = lambda i: -1 if i % 2 else 1
   for fy in (2, 3):
       t = toward(fy)
       Bx2[(3, fy)] += h * (3 * (d1[(3, fy)] - d1[(2, fy)]) + (d1[(3, fy + t)] - d1[(2, fy + t)])) / 8
   for fx in (2, 3):
       t = toward(fx)
       By2[(fx, 3)] += h * (3 * (d1[(fx, 3)] - d1[(fx, 2)]) + (d1[(fx + t, 3)] - d1[(fx + t, 2)])) / 8
   d2 = {c: div(Bx2, By2, c) for c in cells}
   coarse_div = (Bxc[2, 1] - Bxc[1, 1]) + (Byc[1, 2] - Byc[1, 1])
   assert abs(coarse_div - q0) < 1e-12 and all(abs(v - q0) < 1e-12 for v in d2.values())
   from phare_draw import C, strip, text, finish
   from matplotlib.patches import Rectangle
   fig, axs = strip(2, titles=["after stage 1", "after stage 2"], phone=True)
   def draw(ax, d):
       ax.add_patch(Rectangle((0, 0), 2, 2, fill=False, lw=2.5, ec=C["ink"]))
       ax.plot([1, 1], [0, 2], color=C["grey"], lw=1)
       ax.plot([0, 2], [1, 1], color=C["grey"], lw=1)
       for (cx, cy), v in d.items():
           ok = abs(v - q0) < 1e-9
           text(ax, cx - 2 + 0.5, cy - 2 + 0.5, f"div B\n{v:.3f}", color=C["ok"] if ok else C["bad"])
   draw(axs[0], d1); draw(axs[1], d2)
   text(axs[1], 1, 2.3, f"coarse cell: div B = {coarse_div:.3f}", color=C["ink"], fs=10)
   finish(fig, axs, max_h=3.0)

The coarse cell of the figure has a non-zero divergence on purpose, 0.3 (a field made of a
potential plus a term that adds 0.3 to every cell). After stage 1 its four fine cells
disagree (0.349, 0.254, 0.346, 0.251: their mean is 0.3). After stage 2 they all have
exactly 0.3: the refinement *copies* the divergence of the coarse cell, it does not
remove it. For a coarse cell with divergence zero, the four fine cells get zero.

At a **regrid**, the new level is made by copying B from the old fine level where the
two overlap and refining it from the coarse level elsewhere; the refinement only writes
the faces that are still empty
(:src:`src/amr/data/field/refine/composite_field_refiner.hpp#L171-L177`). The copied
regions are made of whole coarse cells (fine boxes are, with ratio 2), and, in the
temporary data that SAMRAI builds for the interpolation, the B copies are rounded out to
whole coarse cells so that neighbouring pieces agree on shared faces
(:src:`src/amr/data/field/refine/coarse_cell_round_out.hpp#L108-L125`). A refinement that
would cut a coarse cell in half stops the run with an error instead
(:src:`src/amr/data/field/refine/magnetic_refine_patch_strategy.hpp#L180-L198`). Where a
refined cell touches copied data, the copied faces are the mean of what the coarse level
holds, because synchronization (below) keeps them so: the junction is consistent.

Ghost cells during the sub-steps
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A fine level needs B in its level ghost cells during every sub-step. They are filled in
this order: the coarse B is first interpolated in time, between its values at the start
and at the end of the coarse step, then the result is refined with the two-stage
operator (SAMRAI does the time interpolation on the coarse data before the refinement:
:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L883-L891`,
:src:`src/amr/messengers/mhd_messenger.hpp#L549-L585`). A time mix of two fields with
the same divergence has that divergence, so the refined ghost cells get the divergence of
their coarse parent, with one exception: the ghost cells that touch the fine interior.
There, the face on the border itself belongs to the fine level and is advanced, not
refined, so those cells need not have the parent's divergence. The test of the repository
that measures the divergence on a refined level drops a margin of ghost cells for the same
reason (:src:`tests/simulator/refinement/test_divb_refinement.py#L178-L180`).

None of this changes the divergence of the advanced cells: ghost values enter E, and E
changes B through a curl.

.. _divb-sync:

Synchronization: what the coarse level takes from the fine one
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

B is never averaged from the fine level to the coarse level, in hybrid or in MHD. At the end
of the coarse step the coarse level *redoes* its B update with the electric field of the
fine level (refluxing, :doc:`/dev/refluxing`):

* **Hybrid:** the fine level has summed a quarter of its averaged electric field after each
  of its four sub-steps; the coarse edges get the mean of the two fine edges along them (or
  the value of the fine node on top); the coarse level recomputes
  :math:`\mathbf B=\mathbf B_\text{old}-\Delta t\,\nabla\times\mathbf E`
  (:src:`src/amr/solvers/solver_ppc.hpp#L277-L288`,
  :src:`src/amr/solvers/solver_ppc.hpp#L309-L321`,
  :src:`src/amr/data/field/coarsening/electric_field_coarsener.hpp#L45-L145`).
* **MHD:** the fine level has summed its stage-weighted edge electric field and fluxes; the
  coarse level repeats one Euler update from its saved start state with them, B through
  Faraday with that E
  (:src:`src/amr/solvers/solver_mhd.hpp#L347-L359`,
  :src:`src/amr/messengers/mhd_messenger.hpp#L357-L362`).

Two consequences:

* **Each coarse cell keeps its divergence**, also those just outside the fine level and those
  under it, because the redo is a curl again (the section on one level).
* **The two levels agree on the border.** The circulation of E around a coarse face is the
  sum of the circulations around the fine faces inside it: the fine edges in the middle
  count twice with opposite signs. With the same weights on every fine edge (¼ for four
  equal sub-steps, each sub-step's share of the coarse step in general), the change
  of the magnetic flux through a coarse face equals the sum of the changes through its
  fine faces. As the refinement starts them equal, they stay equal, and the field made of
  both levels has no flux missing at the border.

.. plot::
   :include-source: False
   :alt: A coarse face made of two fine faces, with E_z on its three fine corners; the coarse corners copy the end ones; the middle one cancels in the sum.

   import matplotlib.pyplot as plt
   from phare_draw import C

   dt = 0.1
   a, b, c = 0.30, 0.10, 0.50                 # time-averaged E_z on the fine corners: bottom, middle, top
   dB_f1 = -dt * (b - a) / 0.5                # fine face 1 (length 1/2)
   dB_f2 = -dt * (c - b) / 0.5
   dB_c = -dt * (c - a) / 1.0                 # coarse face (length 1): its corners copy the fine end corners
   flux_f = dB_f1 * 0.5 + dB_f2 * 0.5
   flux_c = dB_c * 1.0
   assert abs(flux_f - flux_c) < 1e-12
   fig, ax = plt.subplots(figsize=(3.4, 3.6))
   ax.plot([0, 0], [0, 1], color=C["L0"], lw=5, solid_capstyle="butt")
   ax.plot([0.35, 0.35], [0, 0.5], color=C["L1"], lw=4, solid_capstyle="butt")
   ax.plot([0.35, 0.35], [0.5, 1.0], color=C["L1"], lw=4, solid_capstyle="butt")
   for y, v, name in [(0, a, "a"), (0.5, b, "b"), (1.0, c, "c")]:
       ax.plot([0.35], [y], "o", color=C["E"], ms=8)
       ax.text(0.45, y, f"E_z = {v:+.2f}", color=C["E"], fontsize=8, va="center")
   for y, v in [(0, a), (1.0, c)]:
       ax.plot([0], [y], "o", color=C["E"], ms=8)
       ax.annotate("", xy=(0.02, y), xytext=(0.32, y), arrowprops=dict(arrowstyle="-|>", color=C["grey"], lw=1.2))
   ax.text(-0.1, 0.5, "coarse face", rotation=90, ha="right", va="center", fontsize=8, color=C["L0"])
   ax.text(0.45, 0.25, f"ΔΦ {dB_f1 * 0.5:+.3f}", fontsize=8, color=C["B"], ha="left", va="center",
           bbox=dict(fc="white", ec="none", pad=1.0))
   ax.text(0.45, 0.75, f"ΔΦ {dB_f2 * 0.5:+.3f}", fontsize=8, color=C["B"], ha="left", va="center",
           bbox=dict(fc="white", ec="none", pad=1.0))
   ax.text(0.9, -0.22, f"fine faces: {flux_f:+.3f}   coarse face: {flux_c:+.3f}", fontsize=8, ha="center", color=C["ok"])
   ax.set_xlim(-0.4, 1.5); ax.set_ylim(-0.35, 1.15); ax.set_axis_off()
   fig.tight_layout()

The figure uses a face of the coarse level (grey) that contains two fine faces (blue). The
coarse corners copy the end values of the fine level, and the middle corner cancels in the
sum of the two fine changes of flux ΔΦ. This holds exactly for a fixed or an adaptive time
step: every fine sub-step is weighted by its share of the coarse step
(:doc:`/dev/refluxing`).

.. _divb-measured:

Measured: two levels, MHD
-------------------------

A documentation run (:doc:`/validation/mhd`) on the Gaussian-potential field of above, with
a static refined box in 2D, a fixed time step and six coarse-step ends, checks the three
statements. The figure shows, at each coarse-step end:

.. plot::
   :include-source: False
   :alt: Against time: the maximum divergence on levels 0 and 1, the difference between a fine cell's divergence and its parent's, and the mismatch of the flux through the coarse faces.

   import os
   import numpy as np
   import matplotlib.pyplot as plt
   from phare_draw import C

   for base in ("../validation/data", "validation/data", "doc/source/validation/data"):
       f = os.path.join(base, "mhd_divb_amr.txt")
       if os.path.exists(f):
           break
   d = np.loadtxt(f, comments="#")
   fig, ax = plt.subplots(figsize=(3.4, 3.4))
   ax.semilogy(d[:, 0], d[:, 1], "o-", color=C["L0"], label="max div B, level 0", ms=5)
   ax.semilogy(d[:, 0], d[:, 2], "s--", color=C["L1"], label="max div B, level 1", ms=5)
   ax.semilogy(d[:, 0], d[:, 3], "^-", color=C["L2"], label="|fine cell − parent cell|, div B", ms=5)
   ax.semilogy(d[:, 0], d[:, 4], "v-", color=C["B"], label="flux mismatch, coarse faces", ms=5)
   ax.axhline(3e-8, color=C["B"], lw=1, ls=":")
   ax.text(d[0, 0], 4.5e-8, "resolution of B (single precision)", fontsize=7, color=C["B"])
   ax.axhline(1e-5, color=C["L2"], lw=1, ls=":")
   ax.text(d[0, 0], 1.5e-5, "resolution of div B, fine grid", fontsize=7, color=C["L2"])
   ax.set_xlabel("time (ends of coarse steps)", fontsize=9); ax.set_ylabel("value", fontsize=9)
   ax.set_ylim(1e-9, 1e0)
   ax.legend(fontsize=7, frameon=False, loc="center right", bbox_to_anchor=(1.0, 0.62))
   for s in ("top", "right"):
       ax.spines[s].set_visible(False)
   fig.tight_layout()

* The divergence of each level stays at its initial value, 0.0214 (level 0) and 0.0214
  (level 1): the refined level has the divergence of the coarse level.
* A fine cell differs from its parent coarse cell by at most 5·10\ :sup:`-6`, below the
  resolution of a divergence computed from single-precision output on the fine grid
  (about 10\ :sup:`-5`; the two dotted lines of the figure).
* The flux through each coarse face under the fine level differs from the sum of the fluxes
  through its two fine faces by at most 3·10\ :sup:`-8`, the single-precision resolution of
  B, while the coarse :math:`B_x` changes by 1.6·10\ :sup:`-2` over the run.

These runs show the three statements to the resolution of the output files. They do not
show that the divergence stays at round-off (:math:`10^{-16}`), because the output is in single
precision.

Tests and what is not tested
----------------------------

* A test of the repository checks the divergence on a refined level for the **hybrid** model
  (2D, a Harris current sheet, a static box and tagging): the fine level must be below
  :math:`10^{-13}` and not above ten times the coarse level. It needs diagnostics written in
  double precision; on a default build it is skipped
  (:src:`tests/simulator/refinement/test_divb_refinement.py#L216-L240`).
* For **MHD**, a test checks that B is identical on faces shared by two patches (3 levels, 2D)
  (:src:`tests/simulator/test_mhd_shared_faces.py#L131`).
* No test checks the divergence of B across levels for MHD. The divergence of ghost cells is
  not checked anywhere.

See also: :doc:`amr` (the algorithm), :doc:`mhd_time_integration` (the integrators), the
developer pages :doc:`/dev/field_refinement`, :doc:`/dev/time_interpolation` and
:doc:`/dev/refluxing`.
