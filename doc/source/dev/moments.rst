.. _dev-moments:

===================================
Moments: deposition and sum borders
===================================

.. container:: breadcrumb

   You are here: :doc:`/dev/index` › :doc:`/dev/one_time_step` ›
   :doc:`/dev/ppc_order` › Moments

**In one sentence:** each patch :term:`deposits <deposition>` its particles on the grid:
the density of each :term:`population` on the moment nodes and its flux on the E edges, with
the weights that gather E to the particles; the neighbours' parts are then added near patch
borders, and the totals Ohm's law reads are formed.

The time loop keeps the charge density, the mass density and the ion current on the E
edges up to date. The flux on the moment nodes, the bulk velocity and the momentum tensor
are computed only when an output asks for them: :doc:`/dev/output_moments`.

Where a particle deposits
-------------------------

.. plot::
   :include-source: False
   :caption: 1D, order 1, one particle at x = 2.3. Density, F_E,y and F_E,z share the node
      weights (in 1D, E_y and E_z sit on the nodes): one fused kernel. F_E,x goes to the cell
      centres with the weights that gather E_x to this particle.

   from phare_draw import *
   from matplotlib.patches import Rectangle
   H = 1.6
   fig, (a1, a2) = strip(2, ["density, F_E,y, F_E,z: nodes 2 and 3",
                             "F_E,x: centres 1.5 and 2.5, where E_x lives"], height=1.8)
   for a, nodes, color in ((a1, {2: 0.7, 3: 0.3}, C["rho"]), (a2, {1.5: 0.2, 2.5: 0.8}, C["E"])):
       a.plot([-0.3, 5.3], [0, 0], color=C["grey"], lw=1)
       for n in range(6):
           mark(a, "primal", n, 0, color=C["grey"], size=6)
           if n < 5:
               mark(a, "dual", n + 0.5, 0, color=C["grey"], size=6)
           text(a, n, -0.35, str(n), fs=11)
       for x, w in nodes.items():
           a.add_patch(Rectangle((x - 0.15, 0), 0.3, w * H, fc=color, ec="white", lw=0.5))
           text(a, x, w * H + 0.25, f"{w}", color=color)
       particles(a, [(2.3, -0.75)], color=C["ion"])
       text(a, 2.75, -0.75, "x = 2.3", ha="left", color=C["ion"])
   finish(fig, [a1, a2])

* **Reset** density and E-edge flux ``fluxE`` of each population before a push
  (``resetMoments``, :src:`src/core/numerics/moments/moments.hpp#L14-L24`).
* **Deposit** domain and patch ghost particles with ``depositDensityAndFluxE``: density
  += w·S on the nodes; each component of ``fluxE`` += w·v·S with the start index and
  weights that gather the same E component to the particle
  (:src:`src/core/numerics/interpolator/interpolator.hpp#L667-L718`).
* Components whose E sits on the nodes (:math:`E_z` in 2D, :math:`E_y, E_z` in 1D, none
  in 3D) go through the density's fused kernel, switched on at compile time (``writes``,
  :src:`src/core/numerics/interpolator/interpolator.hpp#L285-L297`); the others through
  ``depositOnQty_``, the transpose of the gather
  (:src:`src/core/numerics/interpolator/interpolator.hpp#L802-L845`).

Deposit, then sum over patch borders
------------------------------------

.. plot::
   :include-source: False
   :caption: 1D, patches A (nodes 0-5) and B (nodes 5-10), 2 ghost nodes each.
      Each patch only holds its own particles' share on the nodes both store;
      the border sum gives every copy the full value.

   from phare_draw import *
   from matplotlib.patches import Rectangle
   a_part = {0: 1, 1: 1, 2: 1, 3: 1.0, 4: 0.75, 5: 0.5, 6: 0.15, 7: 0.0}
   b_part = {3: 0.0, 4: 0.25, 5: 0.5, 6: 0.85, 7: 1, 8: 1, 9: 1, 10: 1}
   H = 2.4
   def axis(a):
       a.plot([-0.5, 10.5], [0, 0], color=C["grey"], lw=1)
       for n in range(11):
           mark(a, "primal", n, 0, color=C["grey"], size=5)
           text(a, n, -0.35, str(n), fs=11)
   def band(a, dom, gh, y, name, color):
       a.add_patch(Rectangle((dom[0] - 0.3, y), dom[1] - dom[0] + 0.6, 0.3, fc=color, alpha=0.5, ec="none"))
       a.add_patch(Rectangle((gh[0] - 0.3, y), gh[1] - gh[0] + 0.6, 0.3, fc="none", ec=color, hatch="///", lw=0.8))
       text(a, dom[0] - 0.4, y + 0.15, name, color=color, bold=True, ha="right")
   def bars(a, d, color, bottom=None):
       for n, v in d.items():
           b = (bottom or {}).get(n, 0)
           a.add_patch(Rectangle((n - 0.3, b * H), 0.6, v * H, fc=color, ec="white", lw=0.5))
   fig, (a1, a2, a3) = strip(3, [
       "A deposits its particles (ghost nodes 6-7 hatched)",
       "B deposits its particles (ghost nodes 3-4 hatched)",
       "Border sum: on nodes 3-7 each patch adds the other's part"], height=2.2)
   for a in (a1, a2, a3):
       axis(a)
   band(a1, (0, 5), (6, 7), -0.9, "A", C["L1"]); bars(a1, a_part, C["L1"])
   band(a2, (5, 10), (3, 4), -0.9, "B", C["L2"]); bars(a2, b_part, C["L2"])
   band(a3, (0, 5), (6, 7), -0.9, "A", C["L1"]); band(a3, (5, 10), (3, 4), -1.3, "B", C["L2"])
   bars(a3, {n: a_part.get(n, 0) for n in range(0, 11)}, C["L1"])
   bars(a3, {n: b_part.get(n, 0) for n in range(0, 11)}, C["L2"], bottom=a_part)
   brace(a3, 2.6, 7.4, 2.75, "A and B both hold these sums", color=C["ink"])
   finish(fig, [a1, a2, a3])

* **Level ghosts**: on a refined level each patch adds the level ghost particles it owns,
  Old × (1−α) and New × α (``depositOwned``,
  :src:`src/core/numerics/moments/moments.hpp#L27-L42`; :doc:`/dev/level_ghost_particles`).
* **Sum** over patch borders, per population, density and ``fluxE``: each patch adds the
  neighbours' values on the overlap of the ghost boxes, in place, in one fixed order
  (``fillFluxBorders``, ``fillDensityBorders``,
  :src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L393-L410`).
* **Totals**: ``computeLoopMoments`` forms the charge density Σ q·n, the mass density
  Σ m·n and the ion current on the E edges ``chargeFluxE`` = Σ q·F\ :sub:`E`
  (:src:`src/core/data/ions/ions.hpp#L166-L191`).
* **Ghosts beyond a refined level**: ``fillIonMomentGhosts`` fills the ghosts of the
  charge density and ``chargeFluxE`` (next section).

This runs after each of the two pushes of a step (:doc:`/dev/ppc_order`) and when a level
is created or regridded (:src:`src/amr/level_initializer/hybrid_level_initializer.hpp#L97-L132`).

Ghosts beyond a refined level
-----------------------------

.. plot::
   :include-source: False
   :caption: 1D, order 1, right edge of a refined level (node 6). Particles exist up to one
      cell past the edge: node 7, which Ohm reads, gets only part of its density. Its value
      and that of node 8 (the field ghost nodes) are replaced by values refined from the
      coarser level, in space and time.

   from phare_draw import *
   from matplotlib.patches import Rectangle
   H = 1.4
   fig, (a1, a2) = strip(2, ["after the deposit and the border sums",
                             "fillIonMomentGhosts: ghost nodes refined from the coarser level"],
                         height=2.0)
   for a in (a1, a2):
       a.axvspan(6, 9.4, color=C["ghost"], lw=0)
       a.plot([1.6, 9.4], [0, 0], color=C["grey"], lw=1)
       a.plot([6, 6], [-0.25, 1.9], color=C["ink"], lw=2)
       for n in range(2, 10):
           text(a, n, -0.35, str(n), fs=11)
       text(a, 3.8, 1.95, "refined level", color=C["L1"])
       text(a, 7.8, 1.95, "outside", color=C["ghost_edge"])
   for n in range(2, 7):
       for a in (a1, a2):
           a.add_patch(Rectangle((n - 0.15, 0), 0.3, H, fc=C["rho"], ec="white"))
   particles(a1, [(2.4, -0.8), (3.7, -0.8), (5.2, -0.8)], color=C["ion"], size=45)
   particles(a1, [(6.3, -0.8), (6.8, -0.8)], color=C["ghost_edge"], size=45)
   text(a1, 8.4, -0.8, "no particle here", color=C["grey"])
   a1.add_patch(Rectangle((7 - 0.15, 0), 0.3, 0.5 * H, fc=C["bad"], ec="white"))
   text(a1, 7, 0.5 * H + 0.3, "part", color=C["bad"])
   for n in (7, 8):              # 2 ghost nodes at order 1
       a2.add_patch(Rectangle((n - 0.15, 0), 0.3, H, fc=C["ghost_edge"], ec="white"))
   finish(fig, [a1, a2])

* The level ghost particles (orange, one cell wide at order 1) complete every domain node,
  the edge node 6 included, but not node 7: the particles of the next cell are not there.
* ``fillIonMomentGhosts`` sets the ghost nodes of the charge density and ``chargeFluxE`` to
  NaN, copies them from a same-level neighbour where there is one, and refines the rest from
  the coarser level; domain nodes are not touched
  (:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L413-L444`,
  :src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L1055-L1087`).
* In time, the fill interpolates between the coarse copies ``Nold_`` / ``FEold_`` saved by
  ``prepareStep`` and the coarse values already at the new time, like E
  (:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L843-L851`).
* It runs right after the totals, in ``moveIons_``
  (:src:`src/amr/solvers/solver_ppc.hpp#L659-L665`). The root level is not touched: its
  ghosts are complete after the border sums.

The tricks
----------

Flux deposited where E is gathered
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** Ohm's law conserves energy only if its ion current is the one whose product
with E is the work E does on the particles; the flux on the moment nodes, averaged to the
E edges, is not.
**Fix.** ``fluxE`` is deposited with the gather's own start index and weights, so
:math:`\sum_p w_p\,\mathbf{v}_p\cdot\mathbf{E}(\mathbf{x}_p) = \sum_\text{edges}\mathbf{F}_E\cdot\mathbf{E}`
(:src:`src/core/numerics/interpolator/interpolator.hpp#L667-L677`). :issue:`1577`
(energy-conserving Ohm); maths: :doc:`/numerics/hybrid_ohm`.

Sum in place, in one order
^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** Added in the order the transfers run, which differs from patch to patch, the
copies of a node with three contributions or more differ in the last bits.
**Fix.** The transfers store each contribution; ``BorderSums::fold`` adds them in the order
of the contributing patch boxes, the same on every patch
(:src:`src/amr/messengers/field_border_sum.hpp#L35-L50`). Overlap = the ghost boxes'
intersection: ``FieldGhostInterpOverlapFillPattern`` for density,
``TensorFieldGhostInterpOverlapFillPattern`` for ``fluxE``
(:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L944-L955`). :pr:`1581`
(same bits on every patch).

Leavers are deposited by the patch they leave
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** A particle changing patch must be counted once.
**Fix.** The patch it left deposits it from its patch ghost array; the receiver
grabs it after the sum (:doc:`/dev/grabbing_particles`).

One patch per particle
^^^^^^^^^^^^^^^^^^^^^^

**Problem.** Overlapping same-level patches would deposit the overlap twice
(density about ×2).
**Fix.** Each cell has one owner patch; the ``+=`` sum is unchanged. :pr:`1334`.

Each level ghost particle is deposited once
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** Every patch whose ghost box covers a level ghost particle holds a copy; the
border sum would add each copy.
**Fix.** Only the first of these patches deposits it (``LevelGhostOwnership``), before the
border sums (:doc:`/dev/level_ghost_particles`). :pr:`1580` (level ghosts deposited once).

Time weights only for Old and New
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Only the Old and New level ghost arrays are deposited with 1−α and α
(``fillIonPopMomentGhosts``,
:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L448-L495`);
every other deposit uses weight 1, including pushable level ghosts that enter
the domain in push 1.

A NaN left stops every rank
^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** The refine operators write only nodes that are still NaN; a ghost node no
copy or refinement reaches would keep its NaN, and Ohm reads it.
**Fix.** After the fill, any NaN left in these ghosts throws on every rank (``mpi::any``,
:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L1071-L1086`). Hybrid runs are
periodic; non-periodic hybrid boundaries would need a boundary condition for these moments.

Go deeper
---------

* :doc:`/dev/level_ghost_particles` · :doc:`/dev/grabbing_particles` ·
  :doc:`/dev/output_moments`

Related
-------

* :doc:`/dev/ppc_order` · :doc:`/dev/electromagnetic` · :doc:`/dev/fill_patterns` ·
  :doc:`/numerics/particles` · :doc:`/numerics/hybrid_ohm`
