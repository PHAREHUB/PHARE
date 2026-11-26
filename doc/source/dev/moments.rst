.. _dev-moments:

===================================
Moments: deposition and sum borders
===================================

.. container:: breadcrumb

   You are here: :doc:`/dev/index` › :doc:`/dev/one_time_step` ›
   :doc:`/dev/ppc_order` › Moments

**In one sentence:** the ion :term:`moments` (density and flux of each
:term:`population`, then the totals) are made by :term:`depositing <deposition>`
each patch's particles on the grid nodes, then completed near patch borders by
adding the neighbours' parts, in one fixed order, so that every copy of a node
in a patch domain holds the same bits.

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

After each push, on every patch of a level:

.. plot::
   :include-source: False
   :caption: The moment pipeline after a push. Level ghost particles are
      deposited before the border sums, so the sums carry them to every patch.

   from phare_draw import *
   fig, (a,) = strip(1, phone=True, height=4)
   boxes(a, ["Reset\nper population",
             "Deposit\ndomain + patch ghosts",
             "Level ghosts\nowned ones, Old and New",
             "Border sums\nin place, one order",
             "Totals\nρ, V, charge flux"],
         arrows=[(0, 1), (1, 2), (2, 3), (3, 4)], w=4.6, h=0.95, gap=0.45, fs=12)
   finish(fig, [a])

* **Reset** density and flux of each population
  (``resetMoments``, :src:`src/core/numerics/moments/moments.hpp#L14-L22`).
* **Deposit** domain and patch ghost particles on :term:`primal` nodes:
  density += w·S, flux += w·v·S, S the particle shape
  (``ParticleToMesh``, :src:`src/core/numerics/interpolator/interpolator.hpp#L290-L351`;
  ``IonUpdater::updateAndDepositAll_``, :src:`src/core/numerics/ion_updater/ion_updater.hpp#L309-L312`).
* **Level ghosts**: near a refined level's border, each patch adds
  Old × (1−α) and New × α of the level ghost particles it *owns*
  (``fillIonPopMomentGhosts``, ``depositOwned``,
  :src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L418-L455`;
  ownership: :doc:`/dev/level_ghost_particles`).
* **Sum** over patch borders, per population, in place: on the nodes shared by
  several ghost boxes, each patch adds the others' deposits to its own, in one
  order (``fillFluxBorders``, ``fillDensityBorders``,
  :src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L390-L405`;
  ``BorderSums``, :src:`src/amr/messengers/field_border_sum.hpp#L35-L50`).
* **Totals**: ρ = Σ q·n, V = Σ m·F / Σ m·n, charge flux Σ q·F
  (``updateIons``, :src:`src/core/numerics/ion_updater/ion_updater.hpp#L111-L117`;
  :src:`src/core/data/ions/ions.hpp#L85-L197`).

The order of these calls after a push is in ``SolverPPC::moveIons_``
(:src:`src/amr/solvers/solver_ppc.hpp#L620-L632`, :doc:`/dev/ppc_order`); at
level initialization it is the same, with Old only
(:src:`src/amr/level_initializer/hybrid_level_initializer.hpp#L97-L123`).

The tricks
----------

Store, then add in one order
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. plot::
   :include-source: False
   :caption: A node (red) shared by three patches. Each holds one part of its
      density: a, b, c. Floating-point addition is not associative: two orders
      give two different last digits. Neither order is more exact; what
      matters is that every copy uses the same one, sorted by patch box.

   from phare_draw import *
   a, b, c = 0.1, 0.2, 0.3
   fig, (p1, p2) = strip(2, [
       "Patches A, B, C meet at one node",
       "Its value, added in three orders"], height=2.6, phone=True)
   grid(p1, (0, 0, 8, 6), color=C["light"])
   patch(p1, (0, 0, 4, 3), color=C["L1"]); text(p1, 2, 1.5, "A", color=C["L1"], bold=True, fs=14)
   patch(p1, (4, 0, 8, 3), color=C["L2"]); text(p1, 6, 1.5, "B", color=C["L2"], bold=True, fs=14)
   patch(p1, (0, 3, 8, 6), color=C["ion"]); text(p1, 4, 4.8, "C", color=C["ion"], bold=True, fs=14)
   mark(p1, "primal", 4, 3, color=C["bad"], size=13)
   text(p1, 3.3, 2.4, "a", color=C["L1"], bold=True); text(p1, 4.7, 2.4, "b", color=C["L2"], bold=True)
   text(p1, 4.0, 3.7, "c", color=C["ion"], bold=True)
   rows = [("patch A", f"(a+b)+c = {(a + b) + c!r}", C["bad"]),
           ("patch C: copies differ", f"(b+c)+a = {(b + c) + a!r}", C["bad"]),
           ("all, sorted A, C, B: same", f"(a+c)+b = {(a + c) + b!r}", C["ok"])]
   for k, (who, val, col) in enumerate(rows):
       y = -1.6 * k
       text(p2, 0, y, who, ha="left", color=col, bold=True)
       text(p2, 0, y - 0.7, val, ha="left", family="monospace")
   p2.plot([0, 8], [-4.2, 0.4], alpha=0)
   finish(fig, [p1, p2], align=False)

**Problem.** Where three patches or more share a node (corners, thin patches),
each patch added the others' values in the order its transactions ran, so the
copies of the node differed in the last bits: one node, several values in the
same run. Adding in place would also send a neighbour values that already
hold its own part.
**Fix.** The transactions of the sum :term:`schedule` *store* the contribution
they carry instead of adding it (``BorderSums::add``, one slot per transaction,
:src:`src/amr/messengers/field_border_sum.hpp#L124-L146`;
``FieldBorderSumTransaction``, :src:`src/amr/messengers/field_border_sum.hpp#L319-L333`).
Once the schedule is built, ``finalize`` adds the patch's own contribution
(over the nodes the others reach) and sorts all of them by source box, then
patch id (``lexicographicLess``,
:src:`src/amr/messengers/field_border_sum.hpp#L152-L186`). After the schedule
has run, ``fold`` copies the patch's own values into its slot, zeroes those
nodes, then adds every contribution, its own included, in that order, in place
(:src:`src/amr/messengers/field_border_sum.hpp#L196-L214`). The shared nodes
are the intersection of the two ghost boxes (``FieldGhostInterpOverlapFillPattern``,
``TensorFieldGhostInterpOverlapFillPattern``,
:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L882-L893`;
:doc:`/dev/fill_patterns`).
:pr:`1581` (same order everywhere).

No maximum on the totals
^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** Copies of a node used to differ, and a maximum over the shared
nodes of ρ and of the charge flux made them equal again. A maximum picks the
complete copy only for a quantity that is never negative: the flux has a sign,
and the maximum could put a wrong value on a domain node (:issue:`1575`).
**Fix.** The copies are now equal by construction (one order, level ghosts
deposited before the sums), and the maximum is gone. Ghost nodes of the
moments that no sum reaches may stay incomplete; the code does not rely on
them (comment at
:src:`src/amr/level_initializer/hybrid_level_initializer.hpp#L117-L118`).
:pr:`1583` (ion border max removed).

Level ghosts deposited before the sums
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** A level ghost particle sits in the ghost box of every patch that
covers it, so several patches hold a copy. Deposited after the sums, each patch
only had its own copies: two patches could disagree on a node both store.
**Fix.** One patch, the *owner*, deposits each level ghost particle, before the
sums; the sum adds that deposit to the other patches
(:doc:`/dev/level_ghost_particles`, "One patch deposits each level ghost
particle"). :pr:`1580` (deposit once, before sums).

Leavers are deposited by the patch they leave
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** A particle changing patch must be counted once.
**Fix.** The patch it left deposits it from its patch ghost array; the receiver
grabs it after the sum (:doc:`/dev/grabbing_particles`).

One patch per particle
^^^^^^^^^^^^^^^^^^^^^^

**Problem.** Overlapping same-level patches would deposit the overlap twice
(density about ×2).
**Fix.** Each cell has one owner patch; the border sum is unchanged. :pr:`1334`.

Time weights only for Old and New
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Only the owned Old and New level ghost particles are deposited with 1−α and α
(``fillIonPopMomentGhosts``,
:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L442-L452`);
every other deposit uses weight 1, including pushable level ghosts that enter
the domain in push 1.

Go deeper
---------

* :doc:`/dev/level_ghost_particles` · :doc:`/dev/grabbing_particles` ·
  :doc:`/dev/patch_layout_independence`

Related
-------

* :doc:`/dev/ppc_order` · :doc:`/dev/fill_patterns` · :doc:`/numerics/particles`
