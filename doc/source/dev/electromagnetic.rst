.. _dev-electromagnetic:

======================================
Faraday, Ampère, Ohm and the electrons
======================================

.. container:: breadcrumb

   You are here: :doc:`/dev/index` › :doc:`/dev/one_time_step` › Faraday, Ampère, Ohm

**In one sentence:** in the hybrid model, four small per-patch operators turn B and the ion
moments into J, the electron flux and pressure term, and E; each
:term:`PPC` stage runs them in that order. Physics: :doc:`/theory/hybridpic`; the discrete
Ohm's law and why it conserves energy: :doc:`/numerics/hybrid_ohm`.

One stage, on the grid
----------------------

.. plot::
   :include-source: False
   :caption: What each operator reads (arrows) and writes (highlighted circle), on the
      hybrid Yee grid: B on faces; E, J and the electron flux G on edges; moments on nodes.

   from phare_draw import *
   from matplotlib.patches import Circle

   T = ["Faraday", "Ampère", "Electrons", "Ohm"]
   fig, ax = strip(4, titles=T, phone=True)
   grid(ax[0], (0, 0, 1, 1))
   for a in ax[1:]:
       grid(a, (0, 0, 2, 2))

   # 1 Faraday: Bz(cell centre) <- Ex (bottom/top edges), Ey (left/right edges)
   a = ax[0]
   for (x, y) in [(0.5, 0), (0.5, 1)]:
       arrow(a, (x - 0.22, y), (x + 0.22, y), color=C["E"], lw=2.2, ms=12)
   for (x, y) in [(0, 0.5), (1, 0.5)]:
       arrow(a, (x, y - 0.22), (x, y + 0.22), color=C["E"], lw=2.2, ms=12)
   a.add_patch(Circle((0.5, 0.5), 0.13, fc=C["hi"], ec=C["B"], lw=2, zorder=4))
   mark(a, "primal", 0.5, 0.5, color=C["B"])
   text(a, 0.5, 0.5 - 0.3, "Bz", color=C["B"], bold=True)
   text(a, 0.5, -0.22, "Ex", color=C["E"]); text(a, 1.22, 0.5, "Ey", color=C["E"])

   # 2 Ampere: Jz(node (1,1)) <- Bx (faces x=1), By (faces y=1)
   a = ax[1]
   for (x, y) in [(1, 0.5), (1, 1.5)]:
       arrow(a, (x - 0.22, y), (x + 0.22, y), color=C["B"], lw=2.2, ms=12)
   for (x, y) in [(0.5, 1), (1.5, 1)]:
       arrow(a, (x, y - 0.22), (x, y + 0.22), color=C["B"], lw=2.2, ms=12)
   a.add_patch(Circle((1, 1), 0.13, fc=C["hi"], ec=C["J"], lw=2, zorder=4))
   mark(a, "primal", 1, 1, color=C["J"])
   text(a, 1.25, 1.25, "Jz", color=C["J"], bold=True)
   text(a, 1.25, 0.5, "Bx", color=C["B"]); text(a, 0.5, 1.25, "By", color=C["B"])

   # 3 Electrons: on the Ex edge (0.5,1), Gx = Jix - Jx (both already there, no averaging);
   #   the pressure term on the same edge from n at the nodes (0,1), (1,1)
   a = ax[2]
   for x in (0, 1):
       mark(a, "primal", x, 1, color=C["rho"])
       arrow(a, (x, 1), (0.5 + (x - 0.5) * 0.35, 1), color=C["grey"], lw=1.4, ms=10)
   a.add_patch(Circle((0.5, 1), 0.13, fc=C["hi"], ec=C["J"], lw=2, zorder=4))
   text(a, 0.5, 1.4, "Gx = Jix − Jx", color=C["J"], bold=True)
   text(a, 0.5, 0.6, "−Te ∂x ln n", color=C["rho"])
   text(a, -0.25, 1, "n", color=C["rho"])
   a.plot([-0.5, 2], [1, 1], alpha=0)  # room for the labels

   # 4 Ohm: Ex at (0.5,1) <- G, n and B around it (the ideal term: next figure)
   a = ax[3]
   for x in (0, 1):
       mark(a, "primal", x, 1, color=C["ion"])
       arrow(a, (x, 1), (0.5 + (x - 0.5) * 0.35, 1), color=C["grey"], lw=1.4, ms=10)
   for y in (0.5, 1.5):
       mark(a, "primal", 0.5, y, color=C["B"])
       arrow(a, (0.5, y), (0.5, 1 + (y - 1) * 0.4), color=C["grey"], lw=1.4, ms=10)
   a.add_patch(Circle((0.5, 1), 0.13, fc=C["hi"], ec=C["E"], lw=2, zorder=4))
   text(a, 0.88, 1.22, "Ex", color=C["E"], bold=True)
   text(a, -0.3, 1, "Gz\nn", color=C["ion"], fs=11)
   text(a, 0.95, 0.45, "Bz Gy/n", color=C["B"])
   a.plot([-0.5, 2], [1, 1], alpha=0)  # room for the labels
   finish(fig, ax)

#. **Faraday**: :math:`B_z \mathrel{-}= \Delta t\,(\partial_x E_y - \partial_y E_x)`, domain values only.
#. **Ampère**: :math:`J_z = \partial_x B_y - \partial_y B_x`, on the domain plus one layer, never exchanged.
#. **Electrons**: electron flux :math:`\mathbf{G} = N\mathbf{V}_e = \mathbf{J}_i - \mathbf{J}`
   on the E edges, :math:`\mathbf{J}_i` the ion current deposited there (:doc:`/dev/moments`);
   pressure term :math:`-\nabla P_e/N = -T_e\nabla \ln N` (isothermal), on the E edges too.
#. **Ohm**: :math:`\mathbf{E} = -\mathbf{G}\times\mathbf{B}/N - \nabla P_e/N + \eta\mathbf{J}
   - \nu\nabla^2\mathbf{J}`; the ideal term is built on the B nodes (next section).

* The four calls are in ``predictor1_``, ``predictor2_``, ``corrector_``
  (:src:`src/amr/solvers/solver_ppc.hpp#L409-L531`); which B, E and moments go in at each
  stage is on :doc:`/numerics/hybrid_time_integration`.
* Derivatives come from ``curl``, ``grad``, ``laplacian``
  (:src:`src/core/numerics/differential_operators.hpp`): lazy views evaluated where read,
  primal → dual and dual → primal.
* Averaging to another position: ``GridLayout::project<Source, Destination>``
  (:doc:`/dev/gridlayout`).
* Ghosts of B and E are filled afterwards by the messengers (:doc:`/dev/messengers`).

The ideal term, built on the B nodes
------------------------------------

.. plot::
   :include-source: False
   :caption: The G_y B_z part of E_x, in 2D. 1: G_y, on the four E_y edges, is averaged to
      each B_z node and multiplied by B_z/N there. 2: the two products are averaged back to
      the E_x edge. The G_x B_z part of E_y uses the transposes of these two averages, in
      reverse order (E_x edge to B_z node, B_z node to E_y edge): the two parts then cancel in
      the work G·E.

   from phare_draw import *
   from matplotlib.patches import Circle

   fig, ax = strip(2, titles=["Gy to each Bz node, times Bz/N", "Back to the Ex edge"],
                   phone=True)
   for a in ax:
       grid(a, (0, 0, 1, 2))
   a = ax[0]
   for y in (0.5, 1.5):
       for x in (0, 1):
           mark(a, "primal", x, y, color=C["E"])
           arrow(a, (x, y), (0.5 + (x - 0.5) * 0.35, y), color=C["grey"], lw=1.4, ms=10)
       a.add_patch(Circle((0.5, y), 0.12, fc=C["hi"], ec=C["B"], lw=2, zorder=4))
       text(a, 0.5, y + 0.27, r"$(B_z/N)\,\bar G_y$", color=C["B"])
   text(a, -0.3, 0.5, "Gy", color=C["E"]); text(a, 1.3, 0.5, "Gy", color=C["E"])
   a.plot([-0.5, 1.5], [0, 0], alpha=0)  # room for the labels
   a = ax[1]
   for y in (0.5, 1.5):
       mark(a, "primal", 0.5, y, color=C["B"])
       arrow(a, (0.5, y), (0.5, 1 + (y - 1) * 0.45), color=C["grey"], lw=1.4, ms=10)
   a.add_patch(Circle((0.5, 1), 0.12, fc=C["hi"], ec=C["E"], lw=2, zorder=4))
   text(a, 0.9, 1.0, "Ex", color=C["E"], bold=True)
   finish(fig, ax)

* Each product :math:`(B_c/N)\,\bar G_q` is formed on the :math:`B_c` node, then averaged
  to the E edge:
  :math:`E_x = -P_{B_z\to E_x}(w_z P_{E_y\to B_z} G_y) + P_{B_y\to E_x}(w_y P_{E_z\to B_y} G_z)`,
  :math:`\mathbf{w} = \mathbf{B}/N`, cyclic for :math:`E_y, E_z`
  (``Ohm::ideal_``, :src:`src/core/numerics/ohm/ohm.hpp#L147-L193`).
* :math:`P_{a\to b}` is ``GridLayout::project``: a 2-point average along each direction where
  a and b differ in centering. It also averages a product computed on the fly, which no field
  stores (``core::project``, :src:`src/core/data/grid/gridlayoutdefs.hpp#L56-L69`). N goes to the
  B node the same way.
* In 2D, :math:`B_y` sits on the :math:`E_x` edge: there the second average is the identity.
* Why: :math:`P_{a\to b}` and :math:`P_{b\to a}` are transposes of each other, so the map
  :math:`\mathbf{G}\to\mathbf{E}` is antisymmetric and :math:`\sum \mathbf{G}\cdot\mathbf{E} = 0`
  over the grid: the ideal term does no work on the electrons (:doc:`/numerics/hybrid_ohm`).

Where each operator is
----------------------

* **Faraday** — ``Faraday::operator()``, domain box of each B component
  (:src:`src/core/numerics/faraday/faraday.hpp#L27-L56`).
* **Ampère** — ``Ampere::operator()`` (:src:`src/core/numerics/ampere/ampere.hpp#L35-L63`);
  hybrid uses ``AmpereBox{GrownPhysical, 1}``, the layer Ohm reads
  (:src:`src/amr/solvers/solver_field_evolvers.hpp#L113-L117`).
* **Electrons** — ``Electrons::update``: ``computeFlux`` (G), ``computePressure``,
  ``computePressureTerm`` (:src:`src/core/data/electrons/electrons.hpp#L308-L317`).
  The electron velocity is not stored: ``computeVelocity`` writes
  :math:`\mathbf{V}_e = \mathbf{G}/N` into a field the caller owns
  (:src:`src/core/data/electrons/electrons.hpp#L95-L110`).
* **Pressure term** — ``IsothermalPressureTerm``, called by the electron pressure closure
  (:src:`src/core/numerics/pressure_term/pressure_term.hpp#L24-L66`).
* **Ohm** — ``Ohm::operator()(n, Ge, pressureTerm, B, J, E)`` refuses a G or pressure term
  that is not E-centred (:src:`src/core/numerics/ohm/ohm.hpp#L58-L77`); one function per
  term, ``ideal_``, ``resistive_``, ``hyperresistive_``
  (:src:`src/core/numerics/ohm/ohm.hpp#L147-L270`). ``OhmLevelTransformer`` passes it the
  ion charge density and ``electrons.flux()``
  (:src:`src/amr/solvers/solver_hybrid_field_evolvers.hpp#L29-L41`).

The tricks
----------

An ideal term that does no work on the electrons
   **Problem.** Averaging :math:`\mathbf{V}_e` and B separately to the E edge and multiplying
   there gives :math:`\sum\mathbf{G}\cdot\mathbf{E}\neq 0`: particle noise heats the plasma.
   **Fix.** Products on the B nodes, averaged back with the transposed average
   (:src:`src/core/numerics/ohm/ohm.hpp#L134-L146`). :issue:`1577` (energy-conserving Ohm).

The ion current E works on
   **Problem.** The ions + B energy closes only if :math:`\mathbf{J}_i\cdot\mathbf{E}` is the
   work E does on the particles. :math:`E_x` is gathered with its own weights (dual along x),
   so the flux on the moment nodes, averaged to the edges, is not that current.
   **Fix.** Each population deposits its flux on the E edges with the gather weights
   (``fluxE``, :doc:`/dev/moments`); :math:`\mathbf{G} = \sum_s q_s\mathbf{F}_{E,s} - \mathbf{J}`,
   charge-weighted, not mass-weighted
   (:src:`src/core/data/electrons/electrons.hpp#L77-L93`). :issue:`1466` (electron velocity
   from charge).

Ohm reads one node beyond the box
   **Problem.** ``ideal_`` averages G and N from one node outside the box E is computed on.
   **Fix.** ``computeFlux`` fills G on the box grown by one, where J is valid
   (``evalOnBiggerBox``, :src:`src/core/data/electrons/electrons.hpp#L84-L92`). N and
   :math:`\mathbf{J}_i` there are complete after the border sums; on a refined level they
   are refined from the coarser level (``fillIonMomentGhosts``, :doc:`/dev/moments`).

A pressure term whose curl is zero
   **Problem.** :math:`-\nabla P_e/N` written as a gradient divided by N has a discrete curl
   that is not zero in 2D: density noise creates B (:issue:`1527`).
   **Fix.** The closure computes :math:`-T_e\nabla\ln N`, the discrete gradient of a node
   quantity, whose discrete curl vanishes; Ohm only adds it
   (:src:`src/core/numerics/pressure_term/pressure_term.hpp#L15-L26`). :pr:`1528` (pressure
   term without curl).

Dropping η or ν at no cost
   **Problem.** Testing ``eta > 0`` / ``nu > 0`` in every cell costs time.
   **Fix.** ``Constexprifier`` turns the two flags into template arguments once per call,
   so the compiler builds one loop per combination (:src:`src/core/numerics/ohm/ohm.hpp#L70-L76`).

Hyper-resistivity that follows the mesh
   **Problem.** Constant ν damps grid noise four times more on each finer level; the
   older law ν|B|/n vanished at X points.
   **Fix.** ``hyper_mode="spatial"``: :math:`\nu\,(|\mathbf{B}|/(n+0.1)+1)/4^{\ell}`
   (:src:`src/core/numerics/ohm/ohm.hpp#L240-L270`). :pr:`911` (safer hyper-resistivity).
   MHD has its own law, :math:`\nu\,\Delta x^2(|\mathbf{B}|/\rho+1)` (:doc:`/dev/mhd_solver`).

Go deeper
---------

* :doc:`/numerics/hybrid_ohm` — the discrete Ohm's law, its energy budget.
* :doc:`/dev/moments` — how the ion current on the E edges is deposited.
* :doc:`/dev/gridlayout` — centerings, ``project``, loop boxes.
* :doc:`/dev/fields` — what ``Field`` and ``VecField`` are.

Related
-------

* :doc:`/dev/one_time_step` (up) · :doc:`/dev/ppc_order` · :doc:`/dev/moments` ·
  :doc:`/dev/messengers` · :doc:`/dev/mhd_solver`
