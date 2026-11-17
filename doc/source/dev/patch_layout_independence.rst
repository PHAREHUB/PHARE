.. _dev-patch-layout-independence:

================================
Same result whatever the patches
================================

.. container:: breadcrumb

   You are here: :doc:`/dev/index` › Same result whatever the patches

**In one sentence:** a node stored by several patches must hold the same value
on all of them, and cutting a level into patches differently must not change
the physics: PHARE adds, deposits and keeps particles so that both hold.

.. plot::
   :include-source: False
   :caption: The same level, as one patch and as three. The red node is stored
      once in the first case, three times in the second. Within a run the
      three copies hold the same bits; between the two layouts the values
      agree to round-off.

   from phare_draw import *
   fig, (p1, p2) = strip(2, ["One patch: one copy",
                             "Three patches: three copies"], height=2.0, phone=True)
   grid(p1, (0, 0, 8, 6), color=C["light"]); patch(p1, (0, 0, 8, 6), color=C["ink"])
   mark(p1, "primal", 4, 3, color=C["bad"], size=13)
   grid(p2, (0, 0, 8, 6), color=C["light"])
   patch(p2, (0, 0, 4, 3), color=C["L1"]); text(p2, 2, 1.5, "A", color=C["L1"], bold=True, fs=14)
   patch(p2, (4, 0, 8, 3), color=C["L2"]); text(p2, 6, 1.5, "B", color=C["L2"], bold=True, fs=14)
   patch(p2, (0, 3, 8, 6), color=C["ion"]); text(p2, 4, 4.8, "C", color=C["ion"], bold=True, fs=14)
   mark(p2, "primal", 4, 3, color=C["bad"], size=13)
   text(p2, 3.2, 2.4, "in A", color=C["L1"], bold=True); text(p2, 4.8, 2.4, "in B", color=C["L2"], bold=True)
   text(p2, 4, 3.65, "in C", color=C["ion"], bold=True)
   finish(fig, [p1, p2])

Why it matters
--------------

E and B at a node are computed from these copies: copies that differ make the
patches drift apart, and the run then depends on how the level was cut. That
cut changes at every regrid and with load balancing. And without this
property, comparing two layouts cannot serve as a test.

What is guaranteed
------------------

* **Within one run:** every copy of a node shared by the *domains* of several
  patches, periodic images included, holds the same bits. Ghost nodes may
  differ.
* **Between two layouts** of the same level: the same number of particles per
  level, and the same moments, E and B to round-off, not bit for bit. This
  needs a seeded particle load (``"init": {"seed": ...}``,
  :doc:`/numerics/particles`): only then are a cell's initial particles
  independent of the patches.

Both are checked by ``test_patch_layout_independence.py``: small patches
against one patch per level, 1D and 2D, interpolation orders 1 to 3, patches
just wider than the ghost width, and a run with tagging and regrids (shared
nodes and level 0 particle count only)
(:src:`tests/simulator/test_patch_layout_independence.py#L1-L18`). Level 0 E
and B, ghosts included, are checked by ``test_hybrid_shared_fields.py``
(:src:`tests/simulator/test_hybrid_shared_fields.py#L1-L10`).

Where copies could differ
-------------------------

.. list-table::
   :header-rows: 1
   :widths: 34 40 26

   * - Where
     - What keeps the copies equal
     - Page
   * - A node shared by three patches or more (corners, thin patches):
       floating-point sums depend on the order
     - The border sums add the contributions in one order, sorted by patch box
     - :doc:`/dev/moments` (Store, then add in one order)
   * - A level ghost particle held by several patches
     - One patch, the owner, deposits it, before the border sums
     - :doc:`/dev/level_ghost_particles` (One patch deposits…)
   * - A level ghost particle entering a neighbour's domain, at a re-entrant
       corner of the level border
     - The patch whose domain it is not erases its copy
     - :doc:`/dev/level_ghost_particles` (Entering a neighbour's domain)
   * - Signed charge flux on shared nodes (the former maximum)
     - No maximum: the copies are equal by construction
     - :doc:`/dev/moments` (No maximum on the totals)
   * - Overlapping patches of one level
     - Each cell has one owner patch, which keeps its particles
     - :doc:`/dev/moments` (One patch per particle)
   * - Level 0 B at initialization: ghosts and periodic boundary faces
     - B ghosts filled, then one maximum over ghost overlaps
     - :doc:`/dev/ghost_cells` (Level 0 B ghosts at initialization)

Related
-------

* :doc:`/dev/moments` · :doc:`/dev/level_ghost_particles` ·
  :doc:`/dev/ghost_cells` · :doc:`/dev/messengers`
