============================================
MHD reconstruction: values at the cell faces
============================================

**Goal:** know what the MHD solver reconstructs at a cell face, what each reconstruction
and limiter does, what it costs in ghost cells, and how to choose one.

**In one sentence:** the solver stores one value per cell, but the flux at a face needs two
values, one from each side; the *reconstruction* builds them, and the better it is, the
sharper the result at a given grid, never beyond second order for the scheme as a whole
(:doc:`mhd_finite_volume`).

.. contents:: On this page
   :local:
   :depth: 1

Why reconstruct
---------------

A cell holds one number. Taken as is, the value is flat over the cell and jumps at every
face, even where the flow is smooth. Every jump costs numerical diffusion. A
:term:`reconstruction` replaces the flat value by a better guess at the face, using the
neighbouring cells. The :doc:`Riemann solver <mhd_riemann_solvers>` then turns the left and
right guess into a flux.

What is reconstructed
---------------------

* **Density, velocity (3 components), pressure.** Each is reconstructed separately, with the
  same method (:src:`src/core/numerics/reconstructions/reconstructor.hpp#L19-L33`,
  ``Reconstructor::reconstruct``). Before that, the stored density, momentum, total energy
  and B are turned into velocity and pressure cell by cell
  (:src:`src/core/numerics/primite_conservative_converter/to_primitive_converter.hpp#L23-L28`,
  ``eosEtotToP``).
* **Magnetic field, component along the face normal.** Not reconstructed: it lives on the
  face, both sides take that value
  (:src:`src/core/numerics/reconstructions/reconstructor.hpp#L89-L90`).
* **Magnetic field, the two other components.** They live on other faces. They are first
  averaged to the cell centre (in 1D they already sit there), then reconstructed like a cell
  value (:src:`src/core/numerics/reconstructions/reconstructor.hpp#L36-L51`, ``center_reconstruct``).
* **Current (Hall runs only).** Same treatment as the transverse B
  (:src:`src/core/numerics/godunov_fluxes/godunov_fluxes.hpp#L116-L130`).

Nothing is projected on characteristic variables: the reconstruction works on the primitive
variables, one at a time.

.. admonition:: Side of a face
   :class: note

   For the face between cell :math:`i` and cell :math:`i+1`, the *left* value is built from the
   cells around :math:`i` and the *right* value from the cells around :math:`i+1`
   (:src:`src/core/numerics/reconstructions/linear.hpp#L19-L29`; see the next figure).

The five methods
----------------

``reconstruction`` takes one of five values. The next figure shows which cells each one
reads for the face drawn in the middle.

.. plot::
   :include-source: False
   :alt: For each reconstruction, the cells used to build the left value (blue) and the right value (orange) at one face.

   import matplotlib.pyplot as plt
   from phare_draw import C

   methods = [("Constant", 0, 0, 1), ("Linear", 1, 1, 2), ("WENO3", 1, 1, 2),
              ("WENOZ", 2, 2, 3), ("MP5", 2, 2, 3)]
   # (name, cells to the left of the face used for L, ... for R, ghost layers)
   fig, ax = plt.subplots(figsize=(3.5, 2.9))
   for row, (name, lL, lR, g) in enumerate(methods):
       y = len(methods) - row
       for i in range(-3, 3):
           ax.add_patch(plt.Rectangle((i, y - 0.3), 1, 0.6, fc="white", ec=C["light"], lw=1))
       # face between cell -1 and cell 0: left state reads cells -1-lL .. lL-1, right state cells -lR .. lR
       for i in range(-1 - lL, lL):
           ax.add_patch(plt.Rectangle((i, y - 0.3), 1, 0.3, fc=C["L1"], ec="white", lw=1))
       for i in range(-lR, lR + 1):
           ax.add_patch(plt.Rectangle((i, y), 1, 0.3, fc=C["L2"], ec="white", lw=1))
       ax.text(-3.2, y, f"{name}\n{g} ghost{'s' if g > 1 else ''}", ha="right", va="center", fontsize=8)
   ax.plot([0, 0], [0.5, len(methods) + 0.6], color=C["ink"], lw=1.5)
   ax.text(0, len(methods) + 0.75, "face", ha="center", fontsize=8)
   ax.text(-3, 0.05, "blue: cells read for the left value", color=C["L1"], fontsize=7.5, ha="left")
   ax.text(-3, -0.4, "orange: cells read for the right value", color=C["L2"], fontsize=7.5, ha="left")
   ax.set_xlim(-5.0, 3.2); ax.set_ylim(-0.7, len(methods) + 1.1)
   ax.set_axis_off()
   fig.subplots_adjust(left=0.02, right=0.98, top=0.98, bottom=0.02)

.. list-table::
   :header-rows: 1
   :widths: 40 30 30

   * - ``reconstruction``
     - Ghost layers read
     - Order
   * - ``Constant``
     - 1
     - 1
   * - ``Linear``
     - 2
     - 2
   * - ``WENO3``
     - 2
     - 3
   * - ``WENOZ``
     - 3
     - 5
   * - ``MP5``
     - 3
     - 5

* ``Constant``: the cell value itself
  (:src:`src/core/numerics/reconstructions/constant.hpp#L18-L22`).
* ``Linear``: cell value :math:`\pm` half a limited slope; the slope is
  :math:`\mathrm{limiter}(u_i-u_{i-1},\,u_{i+1}-u_i)`; uses the ``limiter`` option
  (:src:`src/core/numerics/reconstructions/linear.hpp#L47-L67`).
* ``WENO3``: weighted mix of two 2-cell estimates; the weights come from the squared jumps
  with a small constant :math:`\varepsilon=10^{-6}`, not scaled to the data
  (:src:`src/core/numerics/reconstructions/weno3.hpp#L46-L80`).
* ``WENOZ``: weighted mix of three 3-cell estimates; weights from the Jiang–Shu smoothness
  indicators with the WENO-Z correction (:math:`\varepsilon=10^{-40}`)
  (:src:`src/core/numerics/reconstructions/wenoz.hpp#L58-L107`).
* ``MP5``: a fifth-order value, clipped to a monotonicity-preserving range (parameter 4);
  the value is kept unchanged when :math:`(f-u)(f-f_\mathrm{MP})<0` holds, with no tolerance
  (:src:`src/core/numerics/reconstructions/mp5.hpp#L69-L88`).

:term:`WENO` means *weighted essentially non-oscillatory*: the stencil is cut into smaller stencils,
and a cell on the other side of a jump is given almost no weight. The formulas follow
Jiang and Shu (1996) for the smoothness indicators, Borges et al. (2008) for the WENO-Z
weights and Suresh and Huynh (1997) for the monotonicity-preserving limits.

Limiters (``Linear`` only)
--------------------------

The :term:`slope limiter` keeps the slope of ``Linear`` from creating new extrema. With
:math:`a=u_i-u_{i-1}`, :math:`b=u_{i+1}-u_i` and :math:`r=a/b`, the limited slope is
:math:`\varphi(r)\,b`:

.. plot::
   :include-source: False
   :alt: Limiter curves: MinMod and VanLeer, the slope divided by the right-hand difference as a function of r.

   import numpy as np
   import matplotlib.pyplot as plt
   from phare_draw import C
   import mhd_recon as m

   r = np.linspace(-1, 4, 501)
   fig, ax = plt.subplots(figsize=(3.4, 2.8))
   ax.plot(r, [m.minmod(x, 1.0) for x in r], color=C["L2"], lw=2.2, label="MinMod")
   ax.plot(r, [m.van_leer(x, 1.0) for x in r], color=C["L1"], lw=2.2, label="VanLeer")
   ax.axhline(0, color=C["light"], lw=1); ax.axvline(0, color=C["light"], lw=1)
   ax.set_xlabel("r = (u_i − u_(i−1)) / (u_(i+1) − u_i)", fontsize=9)
   ax.set_ylabel("limited slope / (u_(i+1) − u_i)", fontsize=9)
   ax.legend(fontsize=9, frameon=False, loc="lower right")
   for s in ("top", "right"):
       ax.spines[s].set_visible(False)
   fig.tight_layout()

* **MinMod**: zero when the two differences have opposite signs, otherwise the smaller of
  the two (:src:`src/core/numerics/slope_limiters/min_mod.hpp#L15-L25`).
* **VanLeer**: zero when they have opposite signs, otherwise twice their product over their
  sum (:src:`src/core/numerics/slope_limiters/van_leer.hpp#L10`).
  The curve above is its slope as a function of :math:`r`: :math:`(r+|r|)/(1+|r|)`.

``MinMod`` is the more diffusive of the two. At a smooth maximum both give a zero slope, so
``Linear`` is first order there.

What each method does at a jump
-------------------------------

The next figure applies the coded formulas (ported line by line in
``doc/source/_ext/mhd_recon.py``) to one profile: a smooth bump and a jump. The dots are the
cell values, the short bars the left (blue) and right (orange) value at each face. The
closer the two bars are, the smaller the numerical jump the Riemann solver sees.

.. plot::
   :include-source: False
   :alt: A smooth bump and a jump, with the left and right face values built by each reconstruction.

   import numpy as np
   import matplotlib.pyplot as plt
   from phare_draw import C
   import mhd_recon as m

   n = 30
   x = np.arange(n) + 0.5
   u = 0.2 + 0.6 * np.exp(-((x - 8) / 3.0) ** 2) + np.where(x > 20, 0.5, 0.0)
   fig, axes = plt.subplots(6, 1, figsize=(3.4, 8.4), sharex=True, sharey=True)
   for ax, (name, (L, R, g)) in zip(axes.flat, m.METHODS.items()):
       for j in range(4, n - 4):
           ax.plot([x[j] + 0.5 - 0.22, x[j] + 0.5], [L(u, j)] * 2, color=C["L1"], lw=2.2)
           ax.plot([x[j] + 0.5, x[j] + 0.5 + 0.22], [R(u, j)] * 2, color=C["L2"], lw=2.2)
       ax.plot(x[3:-3], u[3:-3], "o", color=C["ink"], ms=3.5)
       ax.set_title(name, fontsize=10, loc="left")
       ax.set_xlim(3.5, 26.5)
       ax.set_yticks([0.2, 0.7])
       ax.set_xticks([])
       for s in ("top", "right"):
           ax.spines[s].set_visible(False)
   fig.tight_layout()

What the figure shows, for this profile only:

* On the smooth bump, the two bars nearly coincide for ``WENOZ`` and ``MP5``, are close for
  ``Linear`` and ``WENO3``, and are far apart for ``Constant``.
* A true step stays a step: no reconstruction smears it (the left and right value differ
  by about 0.5 for all of them). How much diffusion follows is the business of the
  :doc:`Riemann solver <mhd_riemann_solvers>`; every face value at the jump stays between
  the two plateau values (0.2 and 0.7).
* Near the top of the bump, ``WENOZ`` and ``MP5`` give face values above the highest cell
  value (0.805 against 0.784); the other three do not go above it.

Cost: ghost cells and stages
----------------------------

The number of ghost layers a run uses is the stencil width plus one, rounded up to an even
number: 2 for ``Constant`` and 4 for all the other four reconstructions
(:src:`src/core/utilities/ghost_width_calculator.hpp#L64-L68`). The extra layer is kept for the
current computed from B (Ampère's law; the code's own words are in
:src:`src/core/numerics/godunov_fluxes/godunov_fluxes.hpp#L85-L86`). Ghost cells are exchanged
between patches and refined levels (:doc:`/numerics/amr`): the four wide-stencil methods exchange
the same number of layers; they differ in the arithmetic done per face.

Order of each reconstruction, on a small linear wave
-----------------------------------------------------

The orders in the table above are those of the reconstruction alone, on smooth data. They
show up in the full solver only on a **small-amplitude linear wave**. The convergence test
uses an Alfvén wave of amplitude :math:`10^{-6}` over one period and measures the error of
:math:`B_y` for several grids (:doc:`/validation/mhd`). The figure uses the local order between two
successive grids.

.. plot::
   :include-source: False
   :alt: Error of By after one period on a small Alfven wave as a function of the cell size, for each reconstruction.
   :caption: L1 error of :math:`B_y` for a small Alfvén wave. ``WENOZ`` and ``MP5`` are run for ten periods, the others for one, so the vertical positions of the two groups are not comparable; only the slopes are. The two 5th-order curves nearly coincide.

   import numpy as np
   import matplotlib.pyplot as plt
   from phare_draw import C
   import os, glob

   def load(name):
       for base in ("../validation/data", "validation/data", "doc/source/validation/data"):
           f = os.path.join(base, name)
           if os.path.exists(f):
               break
       d = np.array([[float(v) for v in l.split()[:3]] for l in open(f) if l[0] != "#" and l.split()[0].isdigit()])
       return d[:, 1], d[:, 2]

   series = [("Constant", "mhd_convergence_Constant_None_t1.txt", C["grey"], 1),
             ("Linear + VanLeer", "mhd_convergence_Linear_VanLeer_t1.txt", C["L1"], 2),
             ("WENO3", "mhd_convergence_WENO3_None_t1.txt", C["L2"], 3),
             ("WENOZ (t = 10)", "mhd_convergence_WENOZ_None_t10.txt", C["ok"], 5),
             ("MP5 (t = 10)", "mhd_convergence_MP5_None_t10.txt", C["B"], 5)]
   fig, ax = plt.subplots(figsize=(3.4, 3.9))
   for name, f, col, o in series:
       dx, e = load(f)
       if name == "WENO3":
           dx, e = dx[:-1], e[:-1]   # the 1600-cell point is at the float32 floor
       ax.loglog(dx, e, "o-", color=col, label=name, ms=4, lw=1.8)
   dx0 = 0.02
   for o, y0, txt in [(1, 6e-8, "slope 1"), (2, 2e-9, "slope 2"), (3, 1e-10, "slope 3"), (5, 1e-11, "slope 5")]:
       xs = np.array([0.005, 0.02])
       ax.loglog(xs, y0 * (xs / 0.02) ** o, "--", color=C["light"], lw=1.2)
       ax.text(0.0045, y0 * (0.0045 / 0.02) ** o, txt, fontsize=8, color=C["grey"], ha="right", va="center")
   ax.set_xlabel("cell size Δx", fontsize=10); ax.set_ylabel("L1 error of B_y", fontsize=10)
   ax.legend(fontsize=7.5, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=2)
   for s in ("top", "right"):
       ax.spines[s].set_visible(False)
   fig.tight_layout()

Local orders between successive grids, from the coarsest pair to the finest
(Δt = 5·10\ :sup:`-4`, ``SSPRK4_5``, ``Rusanov``; numbers and scripts on the
:doc:`validation page </validation/mhd>`):

* ``Constant``, 50 … 800 cells, one period: 0.80, 0.90, 0.95, 0.97.
* ``Linear`` + ``VanLeer``, 50 … 1600 cells, one period: 1.84, 1.96, 2.00, 2.02, 2.03.
* ``Linear`` + ``MinMod``, 50 … 1600 cells, one period: 1.70, 1.80, 1.86, 1.89, 1.93.
* ``WENO3``, 50 … 800 cells, one period: 2.88, 2.94, 2.97, 2.99. (The 1600-cell point is at
  the precision floor of the output files and is left out.)
* ``WENOZ``, 16 … 64 cells, ten periods: 4.82, 4.84, 4.83, 4.86.
* ``MP5``, 16 … 64 cells, ten periods: 4.63, 4.73, 4.80, 4.85.

.. note::

   The 5th-order methods are shown on coarse grids and a longer time: on finer grids their
   error falls under the precision of the single-precision output files, about
   :math:`10^{-13}`, and cannot be measured. A larger wave amplitude does not show these
   orders either, because the error is then dominated by the nonlinearity of the problem,
   not by the reconstruction. These are properties of the reconstruction on a small linear
   wave. **They are not the order of the scheme**, which is second order at best
   (:doc:`mhd_finite_volume`).

How to choose
-------------

What the reconstruction buys is *less numerical diffusion for a given grid*, not a higher
order of the scheme.

* **Cost.** ``Constant`` uses 2 ghost layers, the others 4 (the same for the four); wider
  stencils cost more arithmetic per face, not more exchanged ghost cells.
* **Jumps.** ``Linear`` with a limiter and ``MP5`` are built to avoid new extrema next to a
  jump, the WENO methods to give little weight to the cells across it. ``Linear`` with
  ``VanLeer`` is the choice of the shock-tube and Orszag–Tang tests of the repository
  (:src:`tests/functional/mhd_shock/common.py#L77-L80`,
  :src:`tests/functional/mhd_orszagtang/orszag_tang.py#L59-L62`).
* **Smooth flows.** On the bump above, ``WENOZ`` and ``MP5`` give the smallest difference
  between the left and right value.
* **Time integrator.** The recommendation in the literature is to pair ``WENOZ`` with a time
  integrator of order three or more (Wang and Spiteri 2007). PHARE does not enforce it; the
  tests that use ``WENOZ`` use ``TVDRK3`` or ``SSPRK4_5`` (:doc:`mhd_time_integration`).

.. seealso::
   :doc:`mhd_finite_volume` (the whole scheme and the order of the scheme),
   :doc:`mhd_riemann_solvers`, :doc:`mhd_time_integration`.
