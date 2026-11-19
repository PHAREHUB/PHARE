=====================================
Ohm's law that conserves energy
=====================================

**Goal:** know why the hybrid solver's Ohm's law is written the way it is on the grid, and what it
guarantees: in a plasma with no resistivity and cold electrons, the energy E takes from the
magnetic field is the energy the ions gain, up to the time-integration error. The physics of Ohm's law is on
:doc:`/theory/hybridpic`; where it sits in a time step, on :doc:`hybrid_time_integration`.

Two conditions are needed, one on Ohm's ideal term and one on the ion current it is given. An
older form of PHARE's Ohm's law broke both, and a quiet uniform plasma slowly gained energy
from particle noise (:issue:`1577` (Ohm energy leak)). This page derives the two conditions and
shows how the grid operators meet them.

The budget
----------

.. plot::
   :include-source: False
   :caption: The energy that leaves one reservoir flows through E into the other. The massless
      electrons only carry it, and they keep none, because E is perpendicular to their flux G.

   import matplotlib.patches as mp
   from phare_draw import *

   fig, (a,) = strip(1, aspect=False, height=2.2)
   def box(x, w, title, sub, col, round_=0.08):
       a.add_patch(mp.FancyBboxPatch((x, 0.15), w, 1.2, boxstyle=f"round,pad=0.02,rounding_size={round_}",
                                     fc="white", ec=col, lw=2.2))
       text(a, x + w / 2, 1.0, title, color=col, bold=True, fs=12)
       text(a, x + w / 2, 0.55, sub, fs=10)
   box(0.0, 2.4, "Ions", "kinetic energy K", C["ion"])
   box(3.3, 2.4, "Electrons", "E = −G × B / n\nG · E = 0", C["grey"], round_=0.5)
   box(6.6, 2.4, "Magnetic field", "W_B = Σ B² / 2", C["B"])
   arrow(a, (2.45, 0.75), (3.25, 0.75), color=C["grey"], style="<|-|>", ms=10)
   arrow(a, (5.75, 0.75), (6.55, 0.75), color=C["grey"], style="<|-|>", ms=10)
   text(a, 1.2, -0.15, "dK/dt = Σ J_E · E", color=C["ion"], fs=10)
   text(a, 7.8, -0.15, "dW_B/dt = −Σ J · E", color=C["B"], fs=10)
   a.set_xlim(-0.2, 9.2); a.set_ylim(-0.4, 1.5)
   finish(fig, [a], align=False)

Take a periodic plasma with :math:`T_e = 0`, no resistivity and no hyper-resistivity, in units
where :math:`e = \mu_0 = 1`. Ohm's law keeps only its ideal term. With the electron flux

.. math::

   \mathbf{G} \equiv n\mathbf{V}_e = \mathbf{J}_i - \mathbf{J}, \qquad
   \mathbf{E} = -\mathbf{V}_e\times\mathbf{B} = -\frac{\mathbf{G}\times\mathbf{B}}{n},

where :math:`\mathbf{J}_i` is the ion current and :math:`\mathbf{J} = \nabla\times\mathbf{B}`.
:math:`\mathbf{G}` lives where :math:`\mathbf{J}` and :math:`\mathbf{E}` live: on the E nodes.

#. **Ions.** A macro-particle of weight :math:`w_p` feels :math:`q(\mathbf{E} + \mathbf{v}\times\mathbf{B})`;
   the magnetic part does no work, so
   :math:`\frac{dK}{dt} = \sum_p q\, w_p\, \mathbf{v}_p\cdot\mathbf{E}(\mathbf{x}_p)`, with
   :math:`\mathbf{E}(\mathbf{x}_p)` the field gathered to the particle.
#. **Magnetic field.** With Faraday and Ampère,

   .. math::

      \frac{dW_B}{dt} = -\sum \mathbf{B}\cdot(\nabla\times\mathbf{E})
      = -\sum \mathbf{E}\cdot(\nabla\times\mathbf{B}) = -\sum \mathbf{J}\cdot\mathbf{E}.

   Moving the curl from E onto B is exact on a periodic Yee grid: the discrete curl of E is the
   transpose of the discrete curl of B (both are differences between neighbouring positions of
   the staggered grid, taken in opposite directions).
#. **Add them.** If the particle work equals :math:`\sum \mathbf{J}_i\cdot\mathbf{E}` for the
   same :math:`\mathbf{J}_i` that Ohm's law is given,

   .. math::

      \frac{d(K + W_B)}{dt} = \sum (\mathbf{J}_i - \mathbf{J})\cdot\mathbf{E} = \sum \mathbf{G}\cdot\mathbf{E},

   the power E gives the electrons, with a minus sign. Massless cold electrons cannot store
   energy, so it must be zero.

So the discrete scheme conserves :math:`K + W_B` when:

.. list-table::
   :header-rows: 1
   :widths: 30 35 35

   * - condition
     - what it constrains
     - how PHARE meets it
   * - **A** · :math:`\sum \mathbf{G}\cdot\mathbf{E} = 0` for every :math:`\mathbf{G}`
     - how the ideal term is written on the grid
     - `Condition A: an antisymmetric ideal term`_
   * - **B** · particle work :math:`= \sum \mathbf{J}_i\cdot\mathbf{E}`
     - which ion current Ohm's law is given
     - `Condition B: the current the particles feel`_

Both are statements in continuous time. The time step adds an error that shrinks when
:math:`\Delta t` is halved; that is all that remains (`What it gives`_).

A built-in alarm: B starts uniform and its mean never changes, so fluctuations can only raise
:math:`W_B`. In a scheme that conserves energy, the ions must therefore cool. Ions that heat in
a quiet uniform plasma are a sign of a leak.

Condition A: an antisymmetric ideal term
----------------------------------------

**Ohm's law as a matrix.** For given B and n the ideal term is linear in G. List every value of
G on the grid in one long vector (:math:`G_x` on every :math:`E_x` node, then :math:`G_y`, then
:math:`G_z`), the same for E, and write :math:`\mathbf{E} = -M\mathbf{G}`. The entry
:math:`M_{ij}` is the weight of G at slot :math:`j` in E at slot :math:`i`. Then

.. math::

   \sum_i G_i E_i = -\sum_{i,j} G_i M_{ij} G_j
   = -\sum_i M_{ii}\,G_i^2 - \sum_{i<j} (M_{ij} + M_{ji})\, G_i G_j .

This is zero for every G if and only if :math:`M^{\mathsf T} = -M`: the influence of G at node
:math:`j` on E at node :math:`i` is minus the influence of G at node :math:`i` on E at node
:math:`j`. At one point the cross product is antisymmetric. On the Yee grid :math:`G_x` and
:math:`G_y` live on different nodes and never meet, so the cancellation has to happen between
pairs of nodes, and depends on where each factor is evaluated.

.. plot::
   :include-source: False
   :caption: One 2D cell, B along z, only the x–y coupling (the 2-point averages are left out).
      G_x lives on the E_x node (bottom edge), G_y on the E_y node (left edge); they meet only
      through the B_z node at the centre. Left: dividing each flux component by the density at
      its own node gives two different factors, and a work proportional to 1/n_x − 1/n_y.
      Right: B_z divided by the density at the B_z node, one factor for both directions: the
      two work terms cancel for any n.

   from phare_draw import *

   fig, ax = strip(2, ["V_e = G / n on the E nodes: does work", "ω = B_z / n on the B_z node: no work"],
                   height=2.6)
   for a, good in zip(ax, (False, True)):
       grid(a, (0, 0, 1, 1))
       for x, y in [(0, 0), (1, 0), (0, 1), (1, 1)]:
           mark(a, "primal", x, y, color=C["rho"], size=7)
       text(a, -0.12, 1.12, "n", color=C["rho"], fs=10)
       mark(a, "dual", 0.5, 0.5, color=C["B"], size=11); text(a, 0.72, 0.62, "B_z", color=C["B"], bold=True)
       mark(a, "primal", 0.5, 0, color=C["E"]); text(a, 0.5, -0.2, "E_x, G_x", color=C["E"])
       mark(a, "primal", 0, 0.5, color=C["E"]); text(a, -0.3, 0.5, "E_y, G_y", color=C["E"], ha="right")
       arrow(a, (0.06, 0.5), (0.4, 0.5), color=C["grey"], lw=1.5, ms=10)
       arrow(a, (0.5, 0.4), (0.5, 0.06), color=C["grey"], lw=1.5, ms=10)
       if good:
           text(a, 1.15, 0.5, "E_x = −ω G_y\nE_y = +ω G_x\nwork = ω(−G_xG_y + G_yG_x) = 0",
                color=C["ok"], ha="left", fs=10)
       else:
           text(a, 1.15, 0.5, "E_x = −B_z G_y / n_y\nE_y = +B_z G_x / n_x\nwork ∝ 1/n_x − 1/n_y ≠ 0",
                color=C["bad"], ha="left", fs=10)
       a.set_xlim(-1.0, 3.6)
   finish(fig, ax)

**The ideal term.** Every factor of the cross product, B and 1/n, is evaluated on the B node
where the two flux components meet. For :math:`E_x` (the other components by cyclic
permutation):

.. math::

   E_x = -P_{B_z\to E_x}\!\left(\omega_z\, P_{E_y\to B_z} G_y\right)
         + P_{B_y\to E_x}\!\left(\omega_y\, P_{E_z\to B_y} G_z\right),
   \qquad \omega = \frac{\mathbf{B}}{\hat n},

where :math:`P_{a\to b}` is the 2-point average that moves a value from the :math:`a` nodes
onto a :math:`b` node, and :math:`\hat n` is the density averaged onto the B node. It is the
x component of :math:`-\mathbf{G}\times\boldsymbol{\omega}`, each product taken on the B node.
A 2-point average is the transpose of the reverse average, :math:`P_{a\to b}^{\mathsf T} = P_{b\to a}`,
and multiplying by :math:`\omega` on the B nodes is a diagonal matrix, so the x–y block
:math:`M_{xy} = P_{B_z\to E_x} D\, P_{E_y\to B_z}` and :math:`M_{yx} = -P_{B_z\to E_y} D\, P_{E_x\to B_z}`
satisfy :math:`M_{yx} = -M_{xy}^{\mathsf T}`: M is antisymmetric for any :math:`n > 0` and any B.

.. admonition:: The same check by hand, in 2D
   :class: note

   Positions in cell :math:`(i, j)`: n at :math:`(i, j)`; :math:`E_x, G_x` at
   :math:`(i+\tfrac12, j)`; :math:`E_y, G_y` at :math:`(i, j+\tfrac12)`; :math:`B_z` at
   :math:`(i+\tfrac12, j+\tfrac12)`. Take the :math:`E_x` node and the :math:`E_y` node that share
   the :math:`B_z` node :math:`(i+\tfrac12, j+\tfrac12)`.

   - :math:`E_x(i+\tfrac12, j)` averages the two :math:`B_z` nodes above and below; at each,
     :math:`G_y` is averaged from its two neighbours. Keeping only the :math:`G_y(i, j+\tfrac12)` term:

     .. math::

        E_x(i+\tfrac12, j) = -\tfrac14\,\omega_z(i+\tfrac12, j+\tfrac12)\, G_y(i, j+\tfrac12) + \dots

   - :math:`E_y(i, j+\tfrac12)` averages the two :math:`B_z` nodes left and right, with a plus sign.
     Keeping only the :math:`G_x(i+\tfrac12, j)` term:

     .. math::

        E_y(i, j+\tfrac12) = +\tfrac14\,\omega_z(i+\tfrac12, j+\tfrac12)\, G_x(i+\tfrac12, j) + \dots

   - In :math:`\sum \mathbf{G}\cdot\mathbf{E}` these give
     :math:`-\tfrac14\omega_z G_x G_y + \tfrac14\omega_z G_y G_x = 0`, with the **same**
     :math:`\omega_z` because both paths go through one :math:`B_z` node. Every neighbouring
     (:math:`E_x`, :math:`E_y`) pair shares exactly one :math:`B_z` node, so the whole sum
     cancels pair by pair.

   The older form computed :math:`\mathbf{V}_e = (\mathbf{F} - \mathbf{J})/n` on the density nodes,
   then averaged :math:`\mathbf{V}_e` and B separately to the E node and multiplied there
   ("average, then multiply"). The two weights of a pair then use different averages of
   :math:`B_z` and of :math:`1/n`, so they are not opposite, and the electrons do work as soon as
   :math:`B/n` varies, which particle noise always makes it do. The unit test below measures a
   relative work of :math:`2\times10^{-2}` (1D), :math:`9\times10^{-3}` (2D) and
   :math:`7\times10^{-4}` (3D) for that form.

The same structure is the energy-conserving form of the vorticity term of the shallow-water
equations on a staggered grid (Sadourny, *J. Atmos. Sci.* 32, 680–689, 1975), with
:math:`\mathbf{G}` in the role of the mass flux and :math:`B/n` in that of the potential
vorticity.

Code: ``Ohm::ideal_`` (:src:`src/core/numerics/ohm/ohm.hpp#L146-L201`). It reads G and n one
node beyond the box E is computed on.

Condition B: the current the particles feel
-------------------------------------------

Each component :math:`E_c` is read at a particle with weights :math:`S^c_i` on that component's
own nodes (primal or dual in each direction, :doc:`grid`). Swapping the two sums of the particle
work:

.. math::

   \sum_p q\,w_p\,v_{p,c}\,E_c(\mathbf{x}_p)
   = \sum_p\sum_i q\,w_p\,v_{p,c}\,S^c_i(\mathbf{x}_p)\,E_{c,i}
   = \sum_i E_{c,i}\underbrace{\Big[\sum_p q\,w_p\,v_{p,c}\,S^c_i(\mathbf{x}_p)\Big]}_{J_{E,c,i}} .

The bracket is a current on the :math:`E_c` nodes, deposited with the weights that gather
:math:`E_c`: the deposit is the transpose of the gather. Handing Ohm's law this current
:math:`\mathbf{J}_E` meets condition B exactly.

The usual moment flux :math:`\mathbf{F}`, deposited on the density nodes with primal weights and
averaged to the E nodes, is not that current:

.. plot::
   :include-source: False
   :caption: One particle at x = i + 0.3, order 1, 1D: its weights on the E_x nodes (dual,
      half-integer). Green: the weights that gather E_x to it, also used to deposit J_E. Red:
      the flux deposited on the integer nodes, then averaged to the E_x nodes: spread over three
      nodes. Same total, different node by node; they agree only for a particle exactly on an
      integer node.

   import numpy as np
   from phare_draw import *

   d = 0.3
   gather = [0.5 - d, 0.5 + d, 0.0]
   smeared = [(1 - d) / 2, 0.5, d / 2]
   fig, (a,) = strip(1, aspect=False, height=2.4)
   xs = np.array([-0.5, 0.5, 1.5])
   a.bar(xs - 0.12, gather, width=0.22, color=C["ok"], label="gather weights (J_E)")
   a.bar(xs + 0.12, smeared, width=0.22, color=C["bad"], label="moment flux averaged to E_x")
   for x, g, s in zip(xs, gather, smeared):
       text(a, x - 0.12, g + 0.06, f"{g:.2f}", color=C["ok"], fs=10)
       text(a, x + 0.12, s + 0.06, f"{s:.2f}", color=C["bad"], fs=10)
   for x in (0, 1):
       mark(a, "primal", x, -0.12, color=C["rho"], size=7); text(a, x, -0.27, ["i", "i+1"][x], color=C["rho"], fs=10)
   for x, lab in zip(xs, ["i−½", "i+½", "i+3/2"]):
       mark(a, "dual", x, -0.12, color=C["E"], size=7); text(a, x, -0.27, lab, color=C["E"], fs=10)
   particles(a, [(d, -0.12)], color=C["ion"], size=40)
   text(a, d, -0.38, "particle", color=C["ion"], fs=10)
   a.plot([-1, 2], [0, 0], color=C["grey"], lw=0.8)
   legend(a, [("gather weights (J_E)", C["ok"], "primal"), ("moment flux averaged to E_x", C["bad"], "primal")],
          loc="upper right", bbox=(1.0, 1.05))
   a.set_xlim(-1, 2.2); a.set_ylim(-0.5, 1.0)
   finish(fig, [a], align=False)

With :math:`\mathbf{J}_\text{used}` the current Ohm's law is given, the budget splits into

.. math::

   \frac{d(K + W_B)}{dt} = \sum (\mathbf{J}_E - \mathbf{J}_\text{used})\cdot\mathbf{E}
   + \underbrace{\sum (\mathbf{J}_\text{used} - \mathbf{J})\cdot\mathbf{E}}_{\sum \mathbf{G}\cdot\mathbf{E}\ =\ 0\ \text{(A)}} .

The first sum vanishes only if :math:`\mathbf{J}_\text{used} = \mathbf{J}_E`. With the averaged
moment flux it is fed by particle noise at the grid scale, and does not shrink with
:math:`\Delta t` (measured in :issue:`1577` (Ohm energy leak)).

So each population deposits, besides its density, a flux on the E nodes, ``fluxE``, with the
gather weights. The E components that sit on the density nodes (:math:`E_z` in 2D, :math:`E_y`
and :math:`E_z` in 1D) have the density's weights, and are deposited together with it. Code:
``Interpolator::depositDensityAndFluxE`` (:src:`src/core/numerics/interpolator/interpolator.hpp#L667-L718`);
the electrons form :math:`\mathbf{G} = \sum_s q_s\,\mathbf{fluxE}_s - \mathbf{J}` on the E nodes
(:src:`src/core/data/electrons/electrons.hpp#L77-L93`). The moment flux on the density nodes
is still what the ``flux`` and ``bulkVelocity`` diagnostics write; it is computed when they are
written (:doc:`/user_guide/diagnostics`).

What it gives
-------------

2D, 24 × 24 cells, :math:`\Delta x = 0.5`, 50 particles per cell, thermal speed 0.3, uniform
density, :math:`|\mathbf{B}| = 1`, :math:`T_e = 0`, no resistivity, one seed, run with this
version of PHARE (``validation/scripts/hybrid_energy.py``, registered in ``run_doc_runs.py``):

.. plot::
   :include-source: False
   :caption: Changes since t = 0, in % of the initial energy, B normal to the plane, Δt = 0.01.
      The ions (green) lose what B (magenta) gains; their sum (black) stays at the
      time-integration error, more than four orders of magnitude below either.

   import os
   import numpy as np
   import matplotlib.pyplot as plt
   from phare_draw import C

   for base in ("../validation/data", "validation/data", "doc/source/validation/data"):
       if os.path.exists(os.path.join(base, "hybrid_energy_z_o1.txt")):
           break
   fig, axes = plt.subplots(1, 2, figsize=(6.5, 2.9), sharey=True)
   for ax, o in zip(axes, (1, 3)):
       d = np.loadtxt(os.path.join(base, f"hybrid_energy_z_o{o}.txt"), comments="#")
       ax.plot(d[:, 0], d[:, 1], color=C["ion"], lw=2, label="ions ΔK")
       ax.plot(d[:, 0], d[:, 2], color=C["B"], lw=2, label="field ΔW_B")
       ax.plot(d[:, 0], d[:, 3], color=C["ink"], lw=2, label="sum")
       ax.axhline(0, color=C["light"], lw=0.8)
       ax.set_title(f"order {o}", fontsize=11, loc="left")
       ax.set_xlabel("t (1/Ω_ci)", fontsize=10)
       for s in ("top", "right"):
           ax.spines[s].set_visible(False)
   axes[0].set_ylabel("% of E0", fontsize=10)
   h, l = axes[0].get_legend_handles_labels()
   fig.legend(h, l, loc="upper center", ncol=3, fontsize=9, frameon=False)
   fig.tight_layout(rect=(0, 0, 1, 0.9))

.. plot::
   :include-source: False
   :caption: The test that separates a spatial error from a time-step error: the gain of ions + B
      at t = 2 (B diagonal in the plane, order 1) with Δt and Δt/2, same particles. It falls when
      Δt is halved: what is left is time-integration error. A spatial leak would not move.

   import os
   import numpy as np
   import matplotlib.pyplot as plt
   from phare_draw import C

   for base in ("../validation/data", "validation/data", "doc/source/validation/data"):
       if os.path.exists(os.path.join(base, "hybrid_energy_xy_o1_dt0.01.txt")):
           break
   dts = (0.01, 0.005)
   g = [abs(np.loadtxt(os.path.join(base, f"hybrid_energy_xy_o1_dt{dt}.txt"), comments="#")[-1, 3]) for dt in dts]
   fig, ax = plt.subplots(figsize=(3.4, 2.6))
   ax.loglog(dts, g, "o-", color=C["L1"], ms=7, lw=1.8)
   for dt, v in zip(dts, g):
       ax.annotate(f"{v:.1e} %", (dt, v), textcoords="offset points", xytext=(6, 4), fontsize=9)
   ax.set_xlabel("Δt", fontsize=10); ax.set_ylabel("|gain| at t = 2, % of E0", fontsize=10)
   ax.set_xlim(0.004, 0.013)
   for s in ("top", "right"):
       ax.spines[s].set_visible(False)
   fig.tight_layout()

Checked by unit tests on the production code:

* the ideal term does no work on the electrons, :math:`\sum\mathbf{G}\cdot\mathbf{E} < 10^{-13}`
  relative, 1D to 3D, orders 1 to 3, with a density from 0.05 to 20; the older form fails the
  same test (:src:`tests/core/numerics/ohm/test_main.cpp#L746-L775`);
* the E-node deposit is the transpose of the gather, to :math:`10^{-13}`, for particles that
  hit both start-index branches; the averaged moment flux fails it
  (:src:`tests/core/numerics/interpolator/test_main.cpp#L840-L912`).

.. admonition:: What it costs
   :class: note

   The E components that do not sit on the density nodes need their own deposit per particle:
   two of three in 2D, all three in 3D. Runs dominated by the particles take longer than with
   the averaged moment flux.

Related
-------

* :doc:`hybrid_time_integration` — where Ohm's law is called in a step.
* :doc:`particles` — shapes, gather and deposit.
* :doc:`/theory/hybridpic` — Ohm's law and the electron model.
* :doc:`/dev/electromagnetic` — the code of Ohm's law and the electrons.
