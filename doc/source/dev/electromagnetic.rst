.. _dev-electromagnetic:

======================================
Faraday, Ampère, Ohm and the electrons
======================================

.. container:: breadcrumb

   You are here: :doc:`/dev/index` › :doc:`/dev/one_time_step` › Faraday, Ampère, Ohm

**In one sentence:** in the hybrid model, four small per-patch operators turn B and the ion
moments into J, the electron velocity and pressure, and E; each
:term:`PPC` stage runs them in that order. Physics: :doc:`/theory/hybridpic`.

One stage, on the grid
----------------------

.. plot::
   :include-source: False
   :caption: What each operator reads (arrows) and writes (highlighted circle), on the
      hybrid Yee grid: B on faces, E and J on edges, moments and electron quantities on nodes.

   from phare_draw import *

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
   a.add_patch(__import__("matplotlib").patches.Circle((0.5, 0.5), 0.13, fc=C["hi"], ec=C["B"], lw=2, zorder=4))
   mark(a, "primal", 0.5, 0.5, color=C["B"])
   text(a, 0.5, 0.5 - 0.3, "Bz", color=C["B"], bold=True)
   text(a, 0.5, -0.22, "Ex", color=C["E"]); text(a, 1.22, 0.5, "Ey", color=C["E"])

   # 2 Ampere: Jz(node (1,1)) <- Bx (left/right... faces x=1), By (faces y=1)
   a = ax[1]
   for (x, y) in [(1, 0.5), (1, 1.5)]:
       arrow(a, (x - 0.22, y), (x + 0.22, y), color=C["B"], lw=2.2, ms=12)
   for (x, y) in [(0.5, 1), (1.5, 1)]:
       arrow(a, (x, y - 0.22), (x, y + 0.22), color=C["B"], lw=2.2, ms=12)
   a.add_patch(__import__("matplotlib").patches.Circle((1, 1), 0.13, fc=C["hi"], ec=C["J"], lw=2, zorder=4))
   mark(a, "primal", 1, 1, color=C["J"])
   text(a, 1.25, 1.25, "Jz", color=C["J"], bold=True)
   text(a, 1.25, 0.5, "Bx", color=C["B"]); text(a, 0.5, 1.25, "By", color=C["B"])

   # 3 Electrons: node (1,1) gets N, F (same node) and Jx averaged from edges (0.5,1),(1.5,1)
   a = ax[2]
   for x in (0.5, 1.5):
       arrow(a, (x - 0.22, 1), (x + 0.22, 1), color=C["J"], lw=2.2, ms=12)
       arrow(a, (x, 1.1), (1 + (x - 1) * 0.25, 1.1), color=C["grey"], lw=1.4, ms=10)
   text(a, 0.5, 0.78, "Jx", color=C["J"]); text(a, 1.5, 0.78, "Jx", color=C["J"])
   a.add_patch(__import__("matplotlib").patches.Circle((1, 1), 0.13, fc=C["hi"], ec=C["ion"], lw=2, zorder=4))
   mark(a, "primal", 1, 1, color=C["ion"])
   text(a, 1, 1.4, "N, F → Ve", color=C["ion"], bold=True)

   # 4 Ohm: Ex at (0.5,1) <- Ve, Pe, n at nodes (0,1),(1,1); Bz at (0.5,0.5),(0.5,1.5); By same place
   a = ax[3]
   for x in (0, 1):
       mark(a, "primal", x, 1, color=C["ion"])
       arrow(a, (x, 1), (0.5 + (x - 0.5) * 0.35, 1), color=C["grey"], lw=1.4, ms=10)
   for y in (0.5, 1.5):
       mark(a, "primal", 0.5, y, color=C["B"])
       arrow(a, (0.5, y), (0.5, 1 + (y - 1) * 0.4), color=C["grey"], lw=1.4, ms=10)
   a.add_patch(__import__("matplotlib").patches.Circle((0.5, 1), 0.13, fc=C["hi"], ec=C["E"], lw=2, zorder=4))
   text(a, 0.88, 1.22, "Ex", color=C["E"], bold=True)
   text(a, -0.25, 1, "Ve\nPe\nn", color=C["ion"], fs=11)
   text(a, 0.82, 0.45, "Bz", color=C["B"])
   finish(fig, ax)

#. **Faraday**: :math:`B_z \mathrel{-}= \Delta t\,(\partial_x E_y - \partial_y E_x)`, domain values only.
#. **Ampère**: :math:`J_z = \partial_x B_y - \partial_y B_x`, on the domain plus one layer, never exchanged.
#. **Electrons**: :math:`\mathbf{V}_e = (\mathbf{F} - \mathbf{J})/N` on nodes (J averaged there);
   :math:`P_e = N T_e` (isothermal).
#. **Ohm**: :math:`E_x = -(V_y B_z - V_z B_y) - \partial_x P_e/n + \eta J_x - \nu\nabla^2 J_x`;
   :math:`B_y` already sits on the :math:`E_x` edge.

* The four calls are in ``predictor1_``, ``predictor2_``, ``corrector_``
  (:src:`src/amr/solvers/solver_ppc.hpp#L410-L532`); which B, E and moments go in at each
  stage is on :doc:`/numerics/time_integration`.
* Derivatives come from ``curl``, ``grad``, ``laplacian``
  (:src:`src/core/numerics/differential_operators.hpp`): lazy views evaluated where read,
  primal → dual and dual → primal.
* Averaging to another position: ``GridLayout::project<Source, Destination>``
  (:doc:`/dev/gridlayout`).
* Ghosts of B and E are filled afterwards by the messengers (:doc:`/dev/messengers`).

Where each operator is
----------------------

* **Faraday** — ``Faraday::operator()``, domain box of each B component
  (:src:`src/core/numerics/faraday/faraday.hpp#L28-L56`).
* **Ampère** — ``Ampere::operator()`` (:src:`src/core/numerics/ampere/ampere.hpp#L35-L63`);
  hybrid uses ``AmpereBox{GrownPhysical, 1}``, the layer Ohm reads
  (:src:`src/amr/solvers/solver_field_evolvers.hpp#L108-L120`).
* **Electrons** — ``Electrons::update``: density, velocity, pressure
  (:src:`src/core/data/electrons/electrons.hpp#L264-L272`).
* **Ohm** — one function per term, ``ideal_``, ``pressure_``, ``resistive_``,
  ``hyperresistive_`` (:src:`src/core/numerics/ohm/ohm.hpp#L101-L289`).

The tricks
----------

Charge flux, not ion bulk velocity
   **Problem.** The ion bulk velocity is mass-weighted: wrong for several species.
   **Fix.** :math:`\mathbf{V}_e = (\sum_s q_s\mathbf{F}_s - \mathbf{J})/N_e`, with
   :math:`N_e` the ion charge density
   (:src:`src/core/data/electrons/electrons.hpp#L70-L98`). :issue:`1466` (electron velocity from charge).

Dropping η or ν at no cost
   **Problem.** Testing ``eta > 0`` / ``nu > 0`` in every cell costs time.
   **Fix.** ``Constexprifier`` turns the two flags into template arguments once per call,
   so the compiler builds one loop per combination (:src:`src/core/numerics/ohm/ohm.hpp#L56-L62`).

Hyper-resistivity that follows the mesh
   **Problem.** Constant ν damps grid noise four times more on each finer level; the
   older law ν|B|/n vanished at X points.
   **Fix.** ``hyper_mode="spatial"``: :math:`\nu\,(|\mathbf{B}|/(n+0.1)+1)/4^{\ell}`
   (:src:`src/core/numerics/ohm/ohm.hpp#L259-L289`). :pr:`911` (safer hyper-resistivity).
   MHD has its own law, :math:`\nu\,\Delta x^2(|\mathbf{B}|/\rho+1)` (:doc:`/dev/mhd_solver`).

.. admonition:: Known issue
   :class: note

   In 2D the isothermal term :math:`-\nabla P_e/n` can create B from density noise: its
   discrete curl is not exactly zero. See :issue:`1527` (pressure term creates B); fix proposed
   in :pr:`1528`.

Go deeper
---------

* :doc:`/dev/gridlayout` — centerings, ``project``, loop boxes.
* :doc:`/dev/fields` — what ``Field`` and ``VecField`` are.

Related
-------

* :doc:`/dev/one_time_step` (up) · :doc:`/dev/ppc_order` · :doc:`/dev/moments` ·
  :doc:`/dev/messengers` · :doc:`/dev/mhd_solver`
