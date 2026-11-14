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
adding the neighbours' parts.

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

* **Reset** density and flux of each population before a push
  (``resetMoments``, :src:`src/core/numerics/moments/moments.hpp#L15-L22`).
* **Deposit** domain and patch ghost particles on :term:`primal` nodes:
  density += w·S, flux += w·v·S, S the particle shape
  (``ParticleToMesh``, :src:`src/core/numerics/interpolator/interpolator.hpp#L290-L351`).
* **Sum** over patch borders, per population: copy to a scratch field, add the
  neighbours' values on the overlap of the two ghost boxes, copy back
  (:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L388-L434`).
* **Level ghosts**: add Old × (1−α) and New × α near a refined level's border
  (:doc:`/dev/level_ghost_particles`).
* **Totals**: ρ = Σ q·n, V = Σ m·F / Σ m·n, charge flux Σ q·F
  (:src:`src/core/data/ions/ions.hpp#L85-L197`); then, on nodes shared by two
  ghost boxes, each patch keeps the larger value of ρ and of the charge flux
  (``fillIonBorders``, ``SetMax``,
  :src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L437-L446`).

This runs after each of the two pushes of a step (:doc:`/dev/ppc_order`).

The tricks
----------

Sum into a scratch copy
^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** Adding in place, B would add A's value that already holds B's part.
**Fix.** Write into ``sumField_`` / ``sumVec_``, read the neighbours' untouched
arrays, copy back (comment at :src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L394-L397`).
Overlap = the two ghost boxes' intersection: ``FieldGhostInterpOverlapFillPattern``
for density, ``TensorFieldGhostInterpOverlapFillPattern`` for flux
(:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L880-L908`).

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

Time weights only for Old and New
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Only the Old and New level ghost arrays are deposited with 1−α and α
(``fillIonPopMomentGhosts``,
:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L457-L495`);
every other deposit uses weight 1, including pushable level ghosts that enter
the domain in push 1.

Why two copies of a node can differ before the max
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Each patch adds the level ghosts of its own layer after the border sum: near a
level border, A's ghost node can miss particles only B holds. Why the maximum
was chosen is not documented (it came with :pr:`1146`).

Go deeper
---------

* :doc:`/dev/level_ghost_particles` · :doc:`/dev/grabbing_particles`

Related
-------

* :doc:`/dev/ppc_order` · :doc:`/dev/fill_patterns` · :doc:`/numerics/particles`
