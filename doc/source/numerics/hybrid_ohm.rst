=====================================
Ohm's law that conserves energy
=====================================

**Goal:** know why the hybrid solver's Ohm's law is written the way it is on the grid, and what it
guarantees: in a plasma with no resistivity and cold electrons, the energy E takes from the
magnetic field is the energy the ions gain, up to the time-integration error. The physics of
Ohm's law is on :doc:`/theory/hybridpic`; where it sits in a time step, on
:doc:`hybrid_time_integration`.

**In short.** Two separate problems made a quiet uniform plasma slowly gain energy from particle
noise (:issue:`1577` (Ohm energy leak)), and each has its own fix:

.. list-table::
   :header-rows: 1
   :widths: 34 33 33

   * -
     - problem 1: E does work on the electrons
     - problem 2: two different ion currents
   * - what goes wrong
     - the ideal term of Ohm's law, written on the grid, is not perpendicular to the electron
       flux, so the field gives energy to massless electrons that cannot hold it
     - the ion current Ohm's law is given is not the current the particles feel when they read
       E, so the ions' gain is not the field's loss
   * - fix
     - every product of the cross product is taken on the B node where the two components meet
       (`Fix 1: one B node per product`_)
     - the ions deposit their current on the E edges with the weights that read E
       (`Fix 2: deposit the current where E is read`_)

The two problems are independent: a scheme can be free of one and still have the other.

Notation
--------

Units with :math:`e = \mu_0 = 1`, a periodic box, one or more ion populations.

.. list-table::
   :header-rows: 1
   :widths: 22 48 30

   * - symbol
     - what it is
     - where it lives (2D)
   * - :math:`n`
     - ion charge density (= electron density)
     - nodes :math:`(i, j)`
   * - :math:`\mathbf{F}`
     - ion charge flux deposited with the density's weights (the ``flux`` diagnostic)
     - nodes :math:`(i, j)`
   * - :math:`\mathbf{J} = \nabla\times\mathbf{B}`
     - total current, from Ampère
     - E edges
   * - :math:`\mathbf{J}_i`
     - the ion current handed to Ohm's law
     - E edges
   * - :math:`\mathbf{G} = n\mathbf{V}_e = \mathbf{J}_i - \mathbf{J}`
     - electron flux (minus the electron current)
     - E edges
   * - :math:`\mathbf{E}`, :math:`\mathbf{B}`
     - fields
     - :math:`E_x` on x-edges :math:`(i+\tfrac12, j)`, :math:`E_y` on y-edges
       :math:`(i, j+\tfrac12)`, :math:`E_z` on nodes; :math:`B_z` at cell centres
       :math:`(i+\tfrac12, j+\tfrac12)`
   * - :math:`K`, :math:`W_B`
     - ion kinetic energy (bulk and thermal, :math:`\sum_p \tfrac12 w_p m v_p^2`), magnetic
       energy (:math:`\sum B^2/2`)
     -
   * - :math:`P_{a\to b}`
     - the 2-point average that moves a value from the :math:`a` positions to a :math:`b`
       position; :math:`P_{a\to b}^{\mathsf T} = P_{b\to a}`
     -

.. plot::
   :include-source: False
   :caption: Where each quantity sits in one 2D cell. The moments (n, F) are on the nodes
      (cell corners), the current, the electron flux and E on the edges, B_z at the centre.

   from phare_draw import *

   fig, (a,) = strip(1, height=2.6, phone=True)
   grid(a, (0, 0, 1, 1))
   for x, y in [(0, 0), (1, 0), (0, 1), (1, 1)]:
       mark(a, "primal", x, y, color=C["rho"], size=8)
   text(a, -0.08, -0.16, "$n$, $\\mathbf{F}$, $E_z$", color=C["rho"], ha="right")
   for x, y in [(0.5, 0), (0.5, 1)]:
       mark(a, "primal", x, y, color=C["E"], size=8)
   text(a, 0.5, -0.17, "$E_x, J_x, G_x$", color=C["E"])
   for x, y in [(0, 0.5), (1, 0.5)]:
       mark(a, "primal", x, y, color=C["E"], size=8)
   text(a, -0.1, 0.5, "$E_y, J_y, G_y$", color=C["E"], ha="right")
   mark(a, "dual", 0.5, 0.5, color=C["B"], size=11)
   text(a, 0.5, 0.66, "$B_z$", color=C["B"], bold=True)
   text(a, 0, 1.18, "$(i, j{+}1)$", fs=10); text(a, 1, -0.17, "$(i{+}1, j)$", fs=10, ha="left")
   a.set_xlim(-1.3, 1.6)
   finish(fig, [a], align=False)

The energy budget
-----------------

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
   text(a, 1.2, -0.15, "dK/dt = particle work", color=C["ion"], fs=10)
   text(a, 7.8, -0.15, "dW_B/dt = −Σ J · E", color=C["B"], fs=10)
   a.set_xlim(-0.2, 9.2); a.set_ylim(-0.4, 1.5)
   finish(fig, [a], align=False)

Take :math:`T_e = 0`, no resistivity and no hyper-resistivity. Ohm's law keeps only its ideal
term:

.. math::

   \mathbf{E} = -\mathbf{V}_e\times\mathbf{B} = -\frac{\mathbf{G}\times\mathbf{B}}{n}.

The electrons feel two forces, :math:`-n\mathbf{E}` and :math:`-n\mathbf{V}_e\times\mathbf{B}`.
Being massless, they balance them: :math:`-n(\mathbf{E} + \mathbf{V}_e\times\mathbf{B}) = 0`.
The magnetic force never does work, so the balance forces :math:`\mathbf{E}\perp\mathbf{V}_e`:
the electric force does no work on the electrons either.

#. **Ions.** A macro-particle of weight :math:`w_p` feels :math:`q(\mathbf{E} + \mathbf{v}\times\mathbf{B})`;
   the magnetic part does no work, so
   :math:`\frac{dK}{dt} = \sum_p q\, w_p\, \mathbf{v}_p\cdot\mathbf{E}(\mathbf{x}_p)`, with
   :math:`\mathbf{E}(\mathbf{x}_p)` the field read at the particle (the *gather*).
#. **Magnetic field.** With Faraday and Ampère,

   .. math::

      \frac{dW_B}{dt} = -\sum \mathbf{B}\cdot(\nabla\times\mathbf{E})
      = -\sum \mathbf{E}\cdot(\nabla\times\mathbf{B}) = -\sum_\text{edges} \mathbf{J}\cdot\mathbf{E}.

   Moving the curl from E onto B is exact on a periodic Yee grid: the discrete curl of E is the
   transpose of the discrete curl of B (both are differences between neighbouring positions of
   the staggered grid, taken in opposite directions).
#. **Add them**, and add and subtract :math:`\sum \mathbf{J}_i\cdot\mathbf{E}`:

   .. math::

      \frac{d(K + W_B)}{dt} =
      \underbrace{\Big[\sum_p q\,w_p\,\mathbf{v}_p\cdot\mathbf{E}(\mathbf{x}_p) - \sum_\text{edges}\mathbf{J}_i\cdot\mathbf{E}\Big]}_{\text{problem 2: must be } 0}
      + \underbrace{\sum_\text{edges}\mathbf{G}\cdot\mathbf{E}}_{\text{problem 1: must be } 0}.

The last sum is the work of E on the electrons (with a minus sign). Massless cold electrons
cannot store energy, so it must be zero; the first bracket says the particles must gain exactly
what the field loses through the current Ohm's law was given. Both are statements in continuous
time; the time step adds an error that shrinks when :math:`\Delta t` is halved (`What it
gives`_).

A built-in alarm: B starts uniform and its mean never changes, so fluctuations can only raise
:math:`W_B`. In a scheme that conserves energy, the ions must therefore cool. Ions that heat in
a quiet uniform plasma are a sign of a leak.

.. note::

   With :math:`T_e > 0`, :math:`\sum\mathbf{G}\cdot\mathbf{E}` is no longer zero: E then balances
   the pressure force too, and its work on the electrons goes into their pressure. :math:`T_e = 0`
   is the case where the budget must close exactly, which makes it the test of the scheme.

Problem 1: E does work on the electrons
---------------------------------------

**Why a cross product can fail on the grid.** For given B and n the ideal term is linear in G.
List every value of G in one long vector (:math:`G_x` on every x-edge, then :math:`G_y`, then
:math:`G_z`), the same for E, and write :math:`\mathbf{E} = -M\mathbf{G}`. Then

.. math::

   \sum_i G_i E_i = -\sum_{i,j} G_i M_{ij} G_j
   = -\sum_i M_{ii}\,G_i^2 - \sum_{i<j} (M_{ij} + M_{ji})\, G_i G_j .

For this to vanish whatever the particle noise makes G, every pair must cancel:
:math:`M_{ji} = -M_{ij}`, the matrix must be antisymmetric. At one point the cross product is
antisymmetric. On the Yee grid :math:`G_x` and :math:`G_y` sit on different edges, so the
cancellation has to happen between pairs of edges, and it depends on how the averages that bring
each factor to each edge are chosen.

**The older form, "average, then multiply".** It formed
:math:`\mathbf{V}_e = (\mathbf{F} - P\,\mathbf{J})/n` on the nodes, then averaged
:math:`\mathbf{V}_e` and B separately to each E edge and multiplied there. Take the ions at rest
(:math:`\mathbf{F} = 0`), so that :math:`\mathbf{G} = -\mathbf{J}` lives on the edges, and follow
one pair, :math:`G_y(i, j+\tfrac12)` and :math:`G_x(i+\tfrac12, j)`:

.. plot::
   :include-source: False
   :caption: The older form, around node (i, j). Each flux component goes to the node (×½),
      is divided by n there, and goes to the other E edge (×½), where it is multiplied by B_z
      averaged at that edge. The two paths average different pairs of B_z values (red squares),
      so their weights are not opposite.

   from phare_draw import *

   fig, ax = strip(2, ["G_y(i, j+½) → E_x(i+½, j)", "G_x(i+½, j) → E_y(i, j+½)"], height=2.4)
   for k, a in enumerate(ax):
       grid(a, (-1, -1, 1, 1))
       mark(a, "primal", 0, 0, color=C["rho"], size=9)
       text(a, -0.1, -0.16, "n(i, j)", color=C["rho"], ha="right", fs=10)
       for b in [(0.5, 0.5), (0.5, -0.5), (-0.5, 0.5), (-0.5, -0.5)]:
           mark(a, "dual", *b, color=C["light"], size=9)
       if k == 0:
           src, dst, bs = (0, 0.5), (0.5, 0), [(0.5, 0.5), (0.5, -0.5)]
           text(a, -0.08, 0.55, "$G_y$", color=C["E"], ha="right")
           text(a, 0.55, -0.15, "$E_x$", color=C["E"], ha="left")
       else:
           src, dst, bs = (0.5, 0), (0, 0.5), [(0.5, 0.5), (-0.5, 0.5)]
           text(a, 0.55, -0.15, "$G_x$", color=C["E"], ha="left")
           text(a, -0.08, 0.55, "$E_y$", color=C["E"], ha="right")
       mark(a, "primal", *src, color=C["E"], size=9)
       mark(a, "primal", *dst, color=C["E"], size=9)
       arrow(a, src, (0.03 * (src[0] > 0), 0.03 * (src[1] > 0)), color=C["grey"], lw=1.6, ms=10)
       arrow(a, (0, 0), dst, color=C["ink"], lw=1.6, ms=10)
       for b in bs:
           mark(a, "dual", *b, color=C["bad"], size=11)
       text(a, 1.25, 0.25, ["weight  −¼ · (1/n) · ½[B_z(i+½, j+½) + B_z(i+½, j−½)]",
                            "weight  +¼ · (1/n) · ½[B_z(i+½, j+½) + B_z(i−½, j+½)]"][k],
            color=C["bad"], ha="left", fs=10)
       a.set_xlim(-1.3, 5.6)
   finish(fig, ax)

The pair adds :math:`(w_1 + w_2)\,G_xG_y` to :math:`\sum\mathbf{G}\cdot\mathbf{E}`, with

.. math::

   w_1 + w_2 = \frac{1}{8\,n(i,j)}\Big[B_z(i-\tfrac12, j+\tfrac12) - B_z(i+\tfrac12, j-\tfrac12)\Big].

The density cancels inside the pair (both paths divide by :math:`n(i,j)`), but the two averages
of :math:`B_z` differ as soon as :math:`B_z` varies across the node, which particle noise always
makes it do. With the ions moving, :math:`\mathbf{F}` enters the same way. The unit test below
measures a relative work :math:`\sum\mathbf{G}\cdot\mathbf{E}` of :math:`2\times10^{-2}` (1D),
:math:`9\times10^{-3}` (2D) and :math:`7\times10^{-4}` (3D) for that form, with random fields.

Problem 2: two different ion currents
-------------------------------------

Each component :math:`E_c` is read at a particle with weights :math:`S^c_m` on that component's
own positions :math:`m` (primal or dual in each direction, :doc:`grid`). Swapping the two sums
of the particle work:

.. math::

   \sum_p q\,w_p\,v_{p,c}\,E_c(\mathbf{x}_p)
   = \sum_p\sum_m q\,w_p\,v_{p,c}\,S^c_m(\mathbf{x}_p)\,E_{c,m}
   = \sum_m E_{c,m}\underbrace{\Big[\sum_p q\,w_p\,v_{p,c}\,S^c_m(\mathbf{x}_p)\Big]}_{J_{E,c,m}} .

The bracket is the ion current the particles feel: a current on the :math:`E_c` positions,
deposited with the weights that read :math:`E_c`. The first bracket of the budget is then
:math:`\sum (\mathbf{J}_E - \mathbf{J}_i)\cdot\mathbf{E}`: it vanishes only if Ohm's law is given
:math:`\mathbf{J}_i = \mathbf{J}_E`.

The older form gave Ohm's law the moment flux :math:`\mathbf{F}`, deposited on the nodes with
the density's weights. Brought to the E edges by an average, it is not :math:`\mathbf{J}_E`:

.. plot::
   :include-source: False
   :caption: One particle at x = i + 0.3, order 1, 1D: its weights on the E_x positions (dual,
      half-integer). Green: the weights that read E_x at the particle, so the current it
      feels. Red: the flux deposited on the integer nodes, then averaged to the E_x positions:
      spread over three positions. Same total, different place by place; they agree only for a
      particle exactly on a node.

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

Summed over many particles the two currents agree at long wavelengths, but not at the grid
scale, where particle noise lives. The difference
:math:`\sum(\mathbf{J}_E - P\mathbf{F})\cdot\mathbf{E}` is fed by that noise and does not
shrink with :math:`\Delta t`. It is a separate problem from the first one: with an ideal term
that does no work, the ions would still gain or lose energy that the field did not lose or gain.

Fix 1: one B node per product
-----------------------------

Every factor of the cross product, B and :math:`1/n`, is evaluated on the B node where the two
flux components meet. For :math:`E_x` (the other components by cyclic permutation):

.. math::

   E_x = -P_{B_z\to E_x}\!\left(\omega_z\, P_{E_y\to B_z} G_y\right)
         + P_{B_y\to E_x}\!\left(\omega_y\, P_{E_z\to B_y} G_z\right),
   \qquad \boldsymbol{\omega} = \frac{\mathbf{B}}{\hat n},

with :math:`\hat n` the density averaged onto the B node. In words, for the :math:`B_z` part in
2D:

#. average :math:`G_y` from its two y-edges, and :math:`G_x` from its two x-edges, to the cell
   centre: :math:`\hat G_x`, :math:`\hat G_y`;
#. at the centre, form :math:`\hat E_x = -\omega_z\hat G_y` and :math:`\hat E_y = +\omega_z\hat G_x`;
#. average :math:`\hat E_x` from the two centres above and below each x-edge, and
   :math:`\hat E_y` from the two centres left and right of each y-edge.

.. plot::
   :include-source: False
   :caption: The same pair as above. Both paths go through the one B_z node they share (green),
      with the same ω = B_z / n̂ there, and opposite signs: the two weights cancel for any n
      and any B.

   from phare_draw import *

   fig, ax = strip(2, ["G_y(i, j+½) → E_x(i+½, j)", "G_x(i+½, j) → E_y(i, j+½)"], height=2.4)
   for k, a in enumerate(ax):
       grid(a, (-1, -1, 1, 1))
       mark(a, "primal", 0, 0, color=C["light"], size=9)
       for b in [(0.5, -0.5), (-0.5, 0.5), (-0.5, -0.5)]:
           mark(a, "dual", *b, color=C["light"], size=9)
       mark(a, "dual", 0.5, 0.5, color=C["ok"], size=12)
       text(a, 0.62, 0.66, "ω", color=C["ok"], bold=True)
       if k == 0:
           src, dst = (0, 0.5), (0.5, 0)
           text(a, -0.08, 0.55, "$G_y$", color=C["E"], ha="right")
           text(a, 0.55, -0.15, "$E_x$", color=C["E"], ha="left")
       else:
           src, dst = (0.5, 0), (0, 0.5)
           text(a, 0.55, -0.15, "$G_x$", color=C["E"], ha="left")
           text(a, -0.08, 0.55, "$E_y$", color=C["E"], ha="right")
       mark(a, "primal", *src, color=C["E"], size=9)
       mark(a, "primal", *dst, color=C["E"], size=9)
       arrow(a, src, (0.46, 0.46), color=C["grey"], lw=1.6, ms=10)
       arrow(a, (0.46, 0.46), dst, color=C["ink"], lw=1.6, ms=10)
       text(a, 1.25, 0.25, ["weight  −½ · ω(i+½, j+½) · ½ = −¼ ω",
                            "weight  +½ · ω(i+½, j+½) · ½ = +¼ ω"][k],
            color=C["ok"], ha="left", fs=10)
       a.set_xlim(-1.3, 5.6)
   finish(fig, ax)

**Why it is antisymmetric.** Multiplying by :math:`\omega` on the B nodes is a diagonal matrix
:math:`D`, and a 2-point average is the transpose of the reverse average. So the x–y block of M
is :math:`M_{xy} = P_{B_z\to E_x} D\, P_{E_y\to B_z}` and
:math:`M_{yx} = -P_{B_z\to E_y} D\, P_{E_x\to B_z} = -M_{xy}^{\mathsf T}`: M is antisymmetric for
any :math:`n > 0` and any B, and :math:`\sum\mathbf{G}\cdot\mathbf{E} = 0` exactly. The same
structure is the energy-conserving form of the vorticity term of the shallow-water equations on
a staggered grid (Sadourny, *J. Atmos. Sci.* 32, 680–689, 1975), with :math:`\mathbf{G}` in the
role of the mass flux and :math:`B/n` in that of the potential vorticity.

Code: ``Ohm::ideal_`` (:src:`src/core/numerics/ohm/ohm.hpp#L147-L193`). The intermediate
:math:`\hat{\mathbf{E}}` is never stored: each E value sums its two B-node terms directly.

The zero is a sum, not a local property
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

On one edge, :math:`G_xE_x` is not zero, and neither is the sum over the four edges of one
cell. Three components never sit at the same place on the Yee grid, so no point can hold
:math:`\mathbf{G}\cdot\mathbf{E} = 0` the way the continuum does. What is zero locally is the
product at each B node: :math:`\hat{\mathbf{G}}\cdot\hat{\mathbf{E}} = -\omega_z\hat G_x\hat G_y
+ \omega_z\hat G_y\hat G_x = 0` at every cell centre (for the :math:`B_z` part).

.. plot::
   :include-source: False
   :caption: How the sum over edges becomes a sum of zeros. 1: E_x on one edge is the average
      of the Ê_x of the two cells above and below it, so the edge's term G_x E_x splits into two
      halves, one for each cell. 2: one cell collects one half from each of its four edges;
      they add up to Ĝ·Ê at its centre, which is zero.

   from phare_draw import *

   fig, ax = strip(2, ["one edge: G_x E_x ≠ 0, split between two cells",
                       "one cell: its four halves add up to Ĝ · Ê = 0"], height=2.6)
   a = ax[0]
   grid(a, (0, 0, 1, 2))
   for c in [(0.5, 0.5), (0.5, 1.5)]:
       mark(a, "dual", *c, color=C["B"], size=10)
   mark(a, "primal", 0.5, 1, color=C["E"], size=10)
   text(a, -0.1, 1.0, "$G_xE_x$", color=C["E"], ha="right")
   arrow(a, (0.5, 1.04), (0.5, 1.42), color=C["J"], lw=1.6, ms=10)
   arrow(a, (0.5, 0.96), (0.5, 0.58), color=C["J"], lw=1.6, ms=10)
   text(a, 1.15, 1.5, "½ $G_x\\hat E_x$(upper cell)", color=C["J"], ha="left", fs=10)
   text(a, 1.15, 0.5, "½ $G_x\\hat E_x$(lower cell)", color=C["J"], ha="left", fs=10)
   a.set_xlim(-1.4, 4.2)
   a = ax[1]
   grid(a, (0, 0, 1, 1))
   mark(a, "dual", 0.5, 0.5, color=C["ok"], size=12)
   for e, d in [((0.5, 0), (0.5, 0.42)), ((0.5, 1), (0.5, 0.58)), ((0, 0.5), (0.42, 0.5)), ((1, 0.5), (0.58, 0.5))]:
       mark(a, "primal", *e, color=C["E"], size=9)
       arrow(a, e, d, color=C["J"], lw=1.6, ms=10)
   text(a, 1.25, 0.5, "$\\hat E_x\\hat G_x + \\hat E_y\\hat G_y$\n$= \\omega(-\\hat G_y\\hat G_x + \\hat G_x\\hat G_y) = 0$",
        color=C["ok"], ha="left", fs=10)
   a.set_xlim(-1.4, 4.2)
   finish(fig, ax, align=False)

So

.. math::

   \sum_\text{edges}\mathbf{G}\cdot\mathbf{E} = \sum_\text{B nodes}\hat{\mathbf{G}}\cdot\hat{\mathbf{E}} = 0,

where regrouping the edge sum into cells is exactly the transpose of the averages. Edge by edge,
energy moves between the electrons and the field at the grid scale, and the neighbouring cell
gives it back: nothing is created, but nothing is zero on its own either. The same holds for
the curl step of the budget, which is exact only summed over the grid; that one is in every Yee
scheme.

Fix 2: deposit the current where E is read
------------------------------------------

Each population deposits, besides its density, a current on the E edges, ``fluxE``, with the
weights that read E: the deposit is the transpose of the gather, and the total
:math:`\mathbf{J}_i = \sum_s q_s\,\mathbf{fluxE}_s` is the :math:`\mathbf{J}_E` of problem 2. The
electrons then form :math:`\mathbf{G} = \mathbf{J}_i - \mathbf{J}` directly on the edges, where
fix 1 needs it. Both brackets of the budget vanish.

* The E components that sit on the nodes (:math:`E_z` in 2D, :math:`E_y` and :math:`E_z` in 1D)
  have the density's weights, so their current is the moment flux itself, deposited together
  with the density. The others need their own deposit: two of three components in 2D, all three
  in 3D. Code: ``Interpolator::depositDensityAndFluxE``
  (:src:`src/core/numerics/interpolator/interpolator.hpp#L702-L753`).
* The electrons: :math:`\mathbf{G} = \sum_s q_s\,\mathbf{fluxE}_s - \mathbf{J}` on the E edges
  (:src:`src/core/data/electrons/electrons.hpp#L77-L93`).
* The moment flux :math:`\mathbf{F}` on the nodes is no longer deposited at every step. It is
  what the ``flux`` and ``bulkVelocity`` diagnostics write, and it is computed when they are
  written (:doc:`/user_guide/diagnostics`, :doc:`/dev/output_moments`).

One more ghost layer of moments
-------------------------------

Fix 1 makes Ohm's law read the moments one node beyond the patch. Take an x-edge on the patch
border: its :math:`E_x` averages the two B nodes above and below, and the one below is outside
the patch. There, :math:`\hat G_y` needs :math:`G_y` on the two y-edges just outside, and
:math:`\hat n` needs the density on the first ghost row of nodes. On a left or right border the
same holds for :math:`E_y`, with :math:`G_x` on the ghost x-edges. The older form read the
moments on the patch's own nodes only.

.. plot::
   :include-source: False
   :caption: The bottom border of a patch (y = 0 is its first row of nodes; orange: outside).
      Top: to compute E_x on the border edge, fix 1 uses the B node below it, which needs G_y on
      two ghost y-edges and n on two ghost nodes (orange). Bottom: the
      older form read V_e only on the two border nodes (and B_z, whose ghosts are always
      filled).

   from phare_draw import *

   fig, ax = strip(2, ["fix 1: E_x on the border edge", "older form: same edge"], height=2.6)
   for k, a in enumerate(ax):
       a.axhspan(-1.25, 0, color=C["ghost"], lw=0)
       grid(a, (0, -1, 1, 1))
       a.plot([-0.6, 1.6], [0, 0], color=C["ink"], lw=2)
       text(a, -0.65, 0.5, "patch", color=C["L1"], ha="right")
       text(a, -0.65, -0.5, "ghost", color=C["ghost_edge"], ha="right")
       mark(a, "primal", 0.5, 0, color=C["E"], size=11)
       text(a, 0.62, 0.2, "$E_x$", color=C["E"], bold=True, ha="left")
       for b in [(0.5, 0.5), (0.5, -0.5)]:
           mark(a, "dual", *b, color=C["B"], size=10)
       if k == 0:
           for p in [(0, 0), (1, 0)]:
               mark(a, "primal", *p, color=C["rho"], size=8)
           for p in [(0, -1), (1, -1)]:
               mark(a, "primal", *p, color=C["ghost_edge"], size=9)
           for p in [(0, -0.5), (1, -0.5)]:
               mark(a, "primal", *p, color=C["ghost_edge"], size=9)
           text(a, 1.75, -0.5, "$G_y$ on 2 ghost y-edges, $n$ on 2 ghost nodes:\none node beyond the patch",
                color=C["ghost_edge"], ha="left", fs=10)
           text(a, 1.75, 0.5, "the B node below the edge is outside", color=C["B"], ha="left", fs=10)
       else:
           for p in [(0, 0), (1, 0)]:
               mark(a, "primal", *p, color=C["rho"], size=10)
           text(a, 1.75, 0.2, "$V_e$ on the two border nodes only", color=C["rho"], ha="left", fs=10)
       a.set_xlim(-1.6, 6.0)
   finish(fig, ax)

That layer must hold the complete moments, as if the particles of the next patch or level were
there:

* **Inside a level**, the neighbouring patch's particles reach it through the patch ghost
  particles and the border sums, as for any other moment (:doc:`/dev/moments`).
* **At the edge of a refined level** nobody deposits the whole layer: the level ghost particles
  are one cell wide at order 1, so the first ghost node gets only part of its density (half of
  it at order 1). ``fillIonMomentGhosts`` refills the ghost nodes of the charge density and of
  ``chargeFluxE`` from the coarser level, in space and time, like E
  (:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L418-L449`; drawing on
  :doc:`/dev/moments`). Without it, the Hall field on the edge row of a refined level would be
  computed from a partial density and flux.

What it gives
-------------

2D, 24 × 24 cells, :math:`\Delta x = 0.5`, 50 particles per cell, thermal speed 0.3, uniform
density, :math:`|\mathbf{B}| = 1`, :math:`T_e = 0`, no resistivity, one seed, run with this
version of PHARE (``validation/scripts/hybrid_energy.py``, registered in ``run_doc_runs.py``;
the data headers name the commit they were run on):

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

Tests
-----

Unit tests on the production code:

* **No work on the electrons**: random positive n (from 0.05 to 20) and random B and G on a
  periodic patch; :math:`|\sum\mathbf{G}\cdot\mathbf{E}| < 10^{-13}` relative to
  :math:`\sum|G_iE_i|`, 1D to 3D, orders 1 to 3; the older form fails the same test
  (:src:`tests/core/numerics/ohm/test_main.cpp#L746-L775`).
* **The deposit is the transpose of the gather**: random particles, including some in the ghost
  cells, with random E and B; the particle work read by the production gather equals
  :math:`\sum \mathbf{fluxE}\cdot\mathbf{E}` from the production deposit to :math:`10^{-13}`,
  for particles that hit both start-index branches in every direction; the moment flux averaged
  to the E edges fails it (:src:`tests/core/numerics/interpolator/test_main.cpp#L877-L949`).
* **The electrons**: :math:`\mathbf{G} = \sum_s q_s\,\mathbf{fluxE}_s - \mathbf{J}` up to one
  node beyond the patch, and with two ion species the electron velocity is charge weighted
  (``ThatElectronFluxIsChargeFluxEMinusJ``, ``ThatElectronVelocityIsChargeWeighted`` in
  ``tests/core/data/electrons/test_electrons.cpp``).

Simulation tests:

* **Ghosts beyond a refined level** (``tests/simulator/test_level_ghost_moments.py``, orders 1
  to 3): a fine ghost node that sits on a coarse node holds exactly the coarse density and
  ``chargeFluxE``, the first ghost node is close to the average of its two coarse neighbours,
  and the flux on the level border is complete; without ``fillIonMomentGhosts`` the ghost nodes
  hold a partial deposit and the test fails.
* **The energy budget** is this page's own run (`What it gives`_), re-run by
  ``run_doc_runs.py`` when the validation data are regenerated.

Related
-------

* :doc:`hybrid_time_integration` — where Ohm's law is called in a step.
* :doc:`particles` — shapes, gather and deposit.
* :doc:`/theory/hybridpic` — Ohm's law and the electron model.
* :doc:`/dev/electromagnetic` — the code of Ohm's law and the electrons.
* :doc:`/dev/moments` — the E-edge deposit, the border sums and the moment ghosts.
