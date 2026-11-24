==================================
MHD Riemann solvers: the face flux
==================================

**Goal:** know what the three Riemann solvers of the MHD module do with the left and right
value at a face, what "more diffusive" means, and how the same solvers give the electric
field on the cell edges.

**In one sentence:** at each face the solver takes the left and right values
(:doc:`mhd_reconstruction`) and returns the flux; the three solvers differ by how many of
the waves leaving the face they keep, hence by how much they smear a discontinuity.

.. contents:: On this page
   :local:
   :depth: 1

Why a solver
------------

After reconstruction, a face has two values, one on each side. If they differ, which
flux crosses the face? The exact answer is the solution of a :term:`Riemann problem`
(two constant states in contact at :math:`t=0`); an approximate :term:`Riemann solver`
gives a cheap estimate from estimates of the wave speeds.

All three solvers write the flux of density, momentum, magnetic field and total energy
through the face. The fluxes of the ideal MHD equations, built from the left and right
values, are the starting point (:src:`src/core/numerics/MHD_equations/MHD_equations.hpp#L32-L45`
for the x direction).

The wave speed they all use
---------------------------

The fastest signal along a face of normal :math:`n` is the fast magnetosonic wave, of speed

.. math::

   c_f^2 = \tfrac12\left(c_s^2 + c_A^2\right)
   + \tfrac12\sqrt{\left(c_s^2 + c_A^2\right)^2 - 4\,c_s^2\,c_{A,n}^2},

with :math:`c_s^2=\gamma P/\rho` the sound speed, :math:`c_A^2=|\mathbf B|^2/\rho` the Alfvén
speed and :math:`c_{A,n}^2=B_n^2/\rho` the Alfvén speed along the normal
(:src:`src/core/numerics/riemann_solvers/mhd_speeds.hpp#L8-L21`). It is computed on the
left and on the right.

The three solvers
-----------------

.. plot::
   :include-source: False
   :alt: Space-time sketch of the wave fan of each solver: Rusanov has one region between two symmetric waves, HLL one region between two waves, HLLD five waves and four star regions.

   import matplotlib.pyplot as plt
   from phare_draw import C

   fig, axes = plt.subplots(3, 1, figsize=(3.4, 7.0))
   def fan(ax, title, speeds, labels, regions):
       for s in speeds:
           ax.plot([0, s], [0, 1], color=C["ink"], lw=1.6)
       ax.axhline(0, color=C["grey"], lw=1.2)
       ax.plot([0], [0], "o", color=C["ink"], ms=4)
       for (x, txt) in regions:
           ax.text(x, 0.16, txt, ha="center", va="center", fontsize=9, color=C["L1"],
                   bbox=dict(fc="white", ec="none", pad=1.0, alpha=0.9), zorder=5)
       for (s, txt) in labels:
           ax.text(s * 1.03, 1.04, txt, ha="center", va="bottom", fontsize=9)
       ax.text(0, -0.1, "face", ha="center", va="top", fontsize=9, color=C["grey"])
       ax.set_title(title, fontsize=10, loc="left")
       ax.set_xlim(-1.35, 1.35); ax.set_ylim(-0.25, 1.28); ax.set_axis_off()
   fan(axes[0], "Rusanov: two waves, speed ∓S", [-1, 1], [(-1, "−S"), (1, "S")],
       [(-1.1, "left"), (0, "one state"), (1.1, "right")])
   fan(axes[1], "HLL: two waves, S_L and S_R", [-1.0, 0.8], [(-1.0, "S_L"), (0.8, "S_R")],
       [(-1.2, "left"), (-0.1, "one state"), (1.1, "right")])
   fan(axes[2], "HLLD: five waves, four inner states",
       [-1.0, -0.55, 0.05, 0.5, 1.0], [(-1.0, "S_L"), (-0.55, "S_L*"), (0.05, "S_M"), (0.5, "S_R*"), (1.0, "S_R")],
       [(-1.25, "L"), (-0.86, "L*"), (-0.36, "L**"), (0.34, "R**"), (0.78, "R*"), (1.25, "R")])
   fig.tight_layout()

Time runs upward from the face. Every wave is a straight line in this picture, between the
left and right state; the solver sets the flux from the states it keeps. The drawing is
schematic: the speeds are not to scale.

* **Rusanov.** One speed :math:`S`, the largest of :math:`|v_n|+c_f` over the left and right
  states. The flux is the average of the two fluxes, minus a diffusion term proportional
  to the jump:
  :math:`F=\tfrac12(F_L+F_R)-\tfrac12\,S\,(U_R-U_L)`
  (:src:`src/core/numerics/riemann_solvers/rusanov.hpp#L85-L110`, flux
  :src:`src/core/numerics/riemann_solvers/rusanov.hpp#L148-L157`).
* **HLL.** Two speeds, the slowest and the fastest signal:
  :math:`S_L=\min(v_{n,L}-c_{f,L},\,v_{n,R}-c_{f,R})`,
  :math:`S_R=\max(v_{n,L}+c_{f,L},\,v_{n,R}+c_{f,R})`. The flux is :math:`F_L` if
  :math:`S_L>0`, :math:`F_R` if :math:`S_R<0`, else
  :math:`\dfrac{S_R F_L - S_L F_R + S_L S_R\,(U_R-U_L)}{S_R-S_L}`
  (:src:`src/core/numerics/riemann_solvers/hll.hpp#L93-L110`, flux
  :src:`src/core/numerics/riemann_solvers/hll.hpp#L154-L170`). Only one state is kept
  between the two waves.
* **HLLD.** The same :math:`S_L` and :math:`S_R`, plus three inner waves: the contact
  (speed :math:`S_M`) and two Alfvén waves (:math:`S_L^*`, :math:`S_R^*`). It keeps four
  inner states, so it can represent a contact discontinuity and the rotation of the
  transverse field across Alfvén waves. The inner states have the form given by Miyoshi and
  Kusano (2005) (:src:`src/core/numerics/riemann_solvers/hlld.hpp#L90-L135`, states
  :src:`src/core/numerics/riemann_solvers/hlld.hpp#L184-L337`).

.. admonition:: A degenerate case of HLLD
   :class: note

   When the Alfvén waves come too close to the fast waves, that is when
   :math:`S_L^*-S_L < 10^{-4}(S_M-S_L)` or :math:`S_R^*-S_R > -10^{-4}(S_R-S_M)`, the
   code sets :math:`S_L^*=S_R^*=S_M` and takes the HLL average for the magnetic field in
   the inner states, so HLLD then behaves like a three-state solver
   (:src:`src/core/numerics/riemann_solvers/hlld.hpp#L208-L220`).

The contribution of the diffusion term is what "diffusive" means: Rusanov always adds
:math:`\tfrac12 S\,(U_R-U_L)`, with the largest speed, to every variable; HLL uses the
two actual speeds; HLLD in addition keeps the contact and the Alfvén waves, so it smears them less.

Hall and dissipative terms
--------------------------

The Riemann solvers solve ideal MHD only: the left and right fluxes, the wave speeds and the
weights left for the edge field contain no Hall, resistive or hyper-resistive term, and the
whistler speed is not in the fan. The whistler speed is used by the adaptive time step only
(:doc:`mhd_time_integration`).

With ``hall=True``, or when ``eta`` or ``nu`` is not zero, the non-ideal electric field
:math:`\mathbf E_\mathrm{ni}` (Hall, resistive and hyper-resistive parts) is computed on the
cell edges from J on the edges (:doc:`mhd_finite_volume`). After the Riemann solve it is added
to the magnetic flux and, as :math:`\mathbf E_\mathrm{ni}\times\mathbf B`, to the energy flux
of each face, and to the upwind electric field on the edges described below.

The electric field on the cell edges
------------------------------------

The magnetic field is advanced with the electric field on the cell edges, not with the
face flux of B (:doc:`mhd_finite_volume`). The edge field is built from the same solvers, in
an *upwind* way (:term:`constrained transport`): it is a weighted average of the values
carried by the faces that touch the edge, plus a diffusion term proportional to the jump of
B.

.. plot::
   :include-source: False
   :alt: A cell corner with the four faces around it, each contributing to the electric field on the corner.

   import matplotlib.pyplot as plt
   from phare_draw import C

   fig, ax = plt.subplots(figsize=(3.4, 2.8))
   for v in (-1, 1):
       ax.plot([v, v], [-1.3, 1.3], color=C["light"], lw=1.2)
       ax.plot([-1.3, 1.3], [v, v], color=C["light"], lw=1.2)
   ax.plot([0, 0], [-1.3, 1.3], color=C["grey"], lw=1.2, ls=":")
   ax.plot([-1.3, 1.3], [0, 0], color=C["grey"], lw=1.2, ls=":")
   ax.plot([0], [0], "o", color=C["E"], ms=10)
   ax.text(0.08, -0.12, "E_z", color=C["E"], fontsize=10, ha="left", va="top", fontweight="bold")
   for (x, y, t, dx, dy) in [(0, 0.55, "x-face above", 0.12, 0), (0, -0.55, "x-face below", 0.12, 0),
                              (0.55, 0, "y-face right", 0, 0.14), (-0.55, 0, "y-face left", 0, 0.14)]:
       ax.plot([x], [y], "s", mfc="white", mec=C["B"], mew=1.8, ms=8)
       ax.annotate("", xy=(0.04 * (1 if x else 0) * (-1 if x < 0 else 1), 0.04 * (1 if y else 0) * (-1 if y < 0 else 1)),
                   xytext=(x, y), arrowprops=dict(arrowstyle="-|>", color=C["B"], lw=1.6))
       ax.text(x + dx, y + dy, t, fontsize=8, color=C["B"], ha="left" if dx else "center",
               va="center" if dx else "bottom")
   ax.set_xlim(-1.3, 1.5); ax.set_ylim(-1.3, 1.3); ax.set_aspect("equal"); ax.set_axis_off()
   fig.tight_layout()

* **Weights.** Each solver leaves, at every face, a left weight, a right weight and two
  diffusion coefficients. For ``Rusanov`` the weights are :math:`\tfrac12,\tfrac12` and the
  coefficients :math:`S/2`. For ``HLL`` the weights are :math:`S_R'/(S_R'-S_L')` and
  :math:`-S_L'/(S_R'-S_L')` and both coefficients are :math:`-S_R'S_L'/(S_R'-S_L')`, with
  :math:`S_L'=\min(0,S_L)` and :math:`S_R'=\max(0,S_R)`.
  For ``HLLD`` they are the HLL ones in the degenerate case above, and weights built from
  the inner speeds otherwise
  (:src:`src/core/numerics/riemann_solvers/rusanov.hpp#L159-L166`,
  :src:`src/core/numerics/riemann_solvers/hll.hpp#L177-L192`,
  :src:`src/core/numerics/riemann_solvers/hlld.hpp#L614-L654`).
  The upwind formulation is that of the constrained-transport literature; the code comment
  cites Mignone and Del Zanna (2021) for the clipped speeds
  (:src:`src/core/numerics/riemann_solvers/hll.hpp#L179-L180`).
* **Velocity at the face.** ``Rusanov`` uses the plain average of the left and right
  velocity. ``HLL`` and ``HLLD`` use the HLL-weighted average. The code comment states that
  this differs from Mignone and Del Zanna and follows the Idefix code (Lesur et al. 2023)
  (:src:`src/core/numerics/riemann_solvers/hll.hpp#L173-L176`).
* **Dimensions.** In 1D the edge field has the form of the flux of B with the sign changed,
  but with one face velocity for both sides. In 2D, :math:`E_x` and :math:`E_y` come from the
  faces of one direction only; :math:`E_z` in 2D, and all three components in 3D, average the
  weights of the faces that meet at the edge
  (:src:`src/core/numerics/constrained_transport/upwind_constrained_transport.hpp#L58-L88`,
  :src:`src/core/numerics/constrained_transport/upwind_constrained_transport.hpp#L91-L114`,
  :src:`src/core/numerics/constrained_transport/upwind_constrained_transport.hpp#L270-L293`).

What "more diffusive" looks like
--------------------------------

The next figure is a *documentation run*, not part of the test suite: the shock tube of
``tests/functional/mhd_shock`` (800 cells, ``TVDRK2``, time step 0.2, time 80) with
``Constant`` reconstruction and each of the three solvers; ``Linear`` with ``VanLeer`` and
``Rusanov`` is added for comparison. Script, command and data are on the
:doc:`validation page </validation/mhd>`.

.. plot::
   :include-source: False
   :alt: Density and B_y in the shock tube at t = 80 for Rusanov, HLL and HLLD with constant reconstruction, and Linear with VanLeer.

   import os
   import numpy as np
   import matplotlib.pyplot as plt
   from phare_draw import C

   def load(name):
       for base in ("../validation/data", "validation/data", "doc/source/validation/data"):
           f = os.path.join(base, name)
           if os.path.exists(f):
               break
       return np.loadtxt(f, comments="#").T

   runs = [("Constant + Rusanov", "mhd_shock_Constant_None_Rusanov.txt", C["bad"]),
           ("Constant + HLL", "mhd_shock_Constant_None_HLL.txt", C["L2"]),
           ("Constant + HLLD", "mhd_shock_Constant_None_HLLD.txt", C["L1"]),
           ("Linear + VanLeer + Rusanov", "mhd_shock_Linear_VanLeer_Rusanov.txt", C["ink"])]
   fig, axes = plt.subplots(2, 1, figsize=(3.4, 5.4), sharex=True)
   for name, f, col in runs:
       x, rho, by = load(f)
       sel = (x >= 380) & (x <= 520)
       axes[0].plot(x[sel], rho[sel], color=col, lw=1.8, label=name)
       axes[1].plot(x[sel], by[sel], color=col, lw=1.8)
   axes[0].set_ylabel("density", fontsize=9); axes[1].set_ylabel("B_y", fontsize=9)
   axes[1].set_xlabel("x", fontsize=9)
   axes[0].legend(fontsize=7, frameon=False, loc="upper right")
   for a in axes:
       for s in ("top", "right"):
           a.spines[s].set_visible(False)
       a.tick_params(labelsize=8)
   fig.tight_layout()

Two waves cross this window. The jump of the density near :math:`x=452` is the **contact
discontinuity**: the density jumps while :math:`B_y` does not (the :math:`B_y` panel is flat
there, at about −0.54 on both sides). A contact does not steepen by itself, so it shows
the diffusion of the solver most clearly. The wave near :math:`x=505` is a slow shock: both
density and :math:`B_y` jump.

The width of each wave, from 10 to 90 % of its density step (the two plateau values are
taken from the sharpest run and used for all, measured by
``doc/source/validation/scripts/mhd_shock_width.py``):

.. list-table::
   :header-rows: 1
   :widths: 60 20 20

   * - Run
     - Contact
     - Slow shock
   * - ``Constant`` + ``Rusanov``
     - 45
     - 12
   * - ``Constant`` + ``HLL``
     - 35
     - 13
   * - ``Constant`` + ``HLLD``
     - 21
     - 7
   * - ``Linear`` + ``MinMod`` + ``Rusanov``
     - 13
     - 4
   * - ``Linear`` + ``VanLeer`` + ``Rusanov``
     - 7
     - 3

Widths are in cells. This is a measure of diffusion on this test, not of accuracy. It shows
that, with the same reconstruction, the solver that keeps the contact wave (``HLLD``) smears
the contact about half as much as ``Rusanov``, and that the reconstruction matters as much.

How to choose
-------------

* The solver sets the *diffusion*: in the run above ``Rusanov`` smears the contact most,
  ``HLLD`` least.
* ``HLL`` and ``HLLD`` are built in 1D with ``Constant`` reconstruction only in the default
  list of compiled modules (:doc:`mhd_finite_volume`).
* No test of the repository checks the result of ``HLL`` or ``HLLD``: they are only
  compiled (:doc:`/validation/mhd`). Every test deck uses ``Rusanov``.

.. seealso::
   :doc:`mhd_reconstruction`, :doc:`mhd_time_integration`, :doc:`mhd_finite_volume`.
