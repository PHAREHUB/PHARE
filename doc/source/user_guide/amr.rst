==================
AMR and refinement
==================

**Goal:** run a hybrid simulation with a finer mesh in part of the box, first
in a region you choose yourself, then in a region PHARE chooses.

AMR stands for *adaptive mesh refinement*: the mesh is made finer only where
the solution needs it, so you do not pay for fine cells everywhere.

.. contents:: On this page
   :local:
   :depth: 1


The idea
--------

The mesh is split in *levels*. Level 0 covers the whole box with the coarse
cells you set with ``cells`` and ``dl``. A finer level covers only part of the
box. Each of its cells is half the size of a level-0 cell in every direction
(the *refinement ratio* is 2, and you cannot change it). A level can itself
have a finer level on top of it, up to 10 levels in all.

A level is not one block but a set of *patches*: boxes of cells that together
cover the refined region.

Finer cells also need shorter time steps. A finer level takes time steps 4
times smaller than the level just below it (the ratio squared), so it does 4
small steps while the coarser level does one.

You choose where the fine levels go in one of two ways:

* ``refinement="boxes"`` (the default): you give the boxes. They never move.
  This is simple and repeatable, and it is the mode of most tests.
* ``refinement="tagging"``: PHARE looks at the solution as it runs, marks
  (*tags*) the cells where the magnetic field changes quickly, and puts fine
  patches around them. The fine regions follow the physics.


Minimal example: boxes you choose
---------------------------------

A 1D run with two levels. Level 0 has 40 cells of size 0.5. The refined region
is the level-0 cells 12 to 27 (both included). Save as ``amr_boxes.py`` in an
empty folder; the setup is the one of :doc:`/getting_started/first_simulation`
(same ``PYTHONPATH``). The run takes a couple of seconds.

.. code-block:: python

   import numpy as np
   import pyphare.pharein as ph
   from pyphare.simulator.simulator import Simulator
   from pyphare.pharesee.run import Run

   # 1. The simulation: 40 cells of size 0.5 on level 0, and one refined region
   sim = ph.Simulation(
       cells=40,
       dl=0.5,                    # the box is 40 * 0.5 = 20 long
       time_step=0.005,
       time_step_nbr=10,
       refinement_boxes={"L0": {"B0": [(12,), (27,)]}},   # level-0 cells 12 to 27
       diag_options={"format": "phareh5",
                     "options": {"dir": "amr_boxes", "mode": "overwrite"}},
   )

   # 2. The initial plasma: uniform density, a bump in By in the middle of the box
   L = 20.0

   def density(x): return 1.0
   def zero(x):    return 0.0
   def bx(x):      return 1.0
   def by(x):      return 0.5 * np.exp(-((x - L / 2) / 2.0) ** 2)
   def vth(x):     return 0.1

   ph.MaxwellianFluidModel(
       bx=bx, by=by, bz=zero,
       protons={"charge": 1, "density": density,
                "vbulkx": zero, "vbulky": zero, "vbulkz": zero,
                "vthx": vth, "vthy": vth, "vthz": vth,
                "nbr_part_per_cell": 50},
   )
   ph.ElectronModel(closure="isothermal", Te=0.01)
   ph.ElectromagDiagnostics(quantity="B", write_timestamps=[0.0, 0.05])

   # 3. Run
   Simulator(sim).run()

   # 4. Look at the patches that were written
   B = Run("amr_boxes").GetB(0.05)
   for ilvl, level in B.levels().items():
       for patch in level.patches:
           print("level", ilvl, "cells", patch.box.lower[0], "to", patch.box.upper[0])

It prints:

.. code-block:: text

   level 0 cells 0 to 39
   level 1 cells 24 to 55

You will also see a warning that some functions are not periodic: the bump
in By is not exactly zero at the box edges. It is harmless here.

Level 1 has one patch. The box you gave covers 16 level-0 cells, so it is
covered by 32 level-1 cells. Their numbers start at 24 because the level-1
cells have their own numbering, see the picture below.


Minimal example: PHARE chooses the boxes
----------------------------------------

Same script, with three changes to the ``Simulation`` call: no
``refinement_boxes``, ``refinement="tagging"``, and the number of levels,
which you must now give yourself (keep the rest of the script as it is).
We also lower the tagging threshold to 0.05
(see the table), because this small bump of By changes slowly: with the default
0.1, nothing would be tagged and the run would have no level 1.

.. code-block:: python

   sim = ph.Simulation(
       cells=40,
       dl=0.5,
       time_step=0.005,
       time_step_nbr=10,
       refinement="tagging",
       max_nbr_levels=2,
       tagging_threshold=0.05,
       diag_options=...,          # unchanged
   )

When we ran it, it printed:

.. code-block:: text

   level 0 cells 0 to 39
   level 1 cells 16 to 63

PHARE found the bump by itself and refined level-0 cells 8 to 31 around it.
Do not rely on these exact numbers: they depend on the threshold and on the
details of the run.


Which cell numbers do I give?
-----------------------------

The boxes in ``refinement_boxes`` are written in the cell numbers of the level
that is being refined, not of the new level. ``"L0"`` boxes use level-0 cell
numbers, and describe where level 1 will be. ``"L1"`` boxes use level-1 cell
numbers (the fine numbering) and describe where level 2 will be. A box is a
pair of cells, the lowest and the highest, both included. The picture shows the
box ``(3,), (6,)`` on level 0.

.. plot::
   :caption: A box on level 0 (cells 3 to 6) and the level-1 cells it becomes
             (cells 6 to 13). Each level-0 cell becomes two level-1 cells:
             cell i becomes cells 2i and 2i+1. Drawn from the rule, not from
             a PHARE run.

   import matplotlib.pyplot as plt
   from matplotlib.patches import Rectangle

   n0, lo, hi = 9, 3, 6                   # level-0 cells, box = cells lo..hi
   fig, ax = plt.subplots(figsize=(5.4, 2.4))
   for i in range(n0):
       inbox = lo <= i <= hi
       ax.add_patch(Rectangle((i, 1.0), 1, 0.8, fc="#cfe3f5" if inbox else "#eeeeee",
                              ec="0.3"))
       ax.text(i + 0.5, 1.4, str(i), ha="center", va="center", fontsize=13)
   for j in range(2 * lo, 2 * hi + 2):
       ax.add_patch(Rectangle((j / 2, 0.0), 0.5, 0.8, fc="#f6d6a8", ec="0.3"))
       ax.text(j / 2 + 0.25, 0.4, str(j), ha="center", va="center", fontsize=11)
   ax.text(-0.15, 1.4, "level 0", ha="right", va="center", fontsize=13)
   ax.text(-0.15, 0.4, "level 1", ha="right", va="center", fontsize=13)
   ax.text((lo + hi + 1) / 2, 2.05, "box (3,), (6,)", ha="center", fontsize=13)
   ax.annotate("", xy=(lo, 1.95), xytext=(hi + 1, 1.95),
               arrowprops=dict(arrowstyle="<->", color="0.2"))
   ax.set_xlim(-2.0, n0 + 0.1)
   ax.set_ylim(-0.1, 2.4)
   ax.axis("off")
   fig.tight_layout()

In 2D or 3D a box is a pair of tuples, ``[(x0, y0), (x1, y1)]``; the same
rule holds in each direction.

The box must stay inside the coarser level, away from its edge: see
``nesting_buffer`` below.


Options
-------

All are arguments of ``ph.Simulation`` (see :doc:`simulation`). Defaults are
read from ``pyphare/pharein/simulation.py``.

.. list-table::
   :header-rows: 1
   :widths: 22 48 30

   * - Option
     - Meaning
     - Default and checks
   * - ``refinement``
     - ``"boxes"``: fixed boxes you give. ``"tagging"``: PHARE chooses the
       boxes during the run. Only these two values work. Python does not
       check the spelling: any other value passes the checks, but no cell is
       ever tagged (see "Common mistakes").
     - ``"boxes"``
   * - ``refinement_boxes``
     - Dictionary, level by level. Keys ``"L0"``, ``"L1"``, ... (an integer key
       such as ``0`` also works). Each value is a dictionary from a name you
       choose (``"B0"``, ...) to a pair ``[(lower), (upper)]`` of cell
       numbers, both included. Each level you list creates the level above
       it. Used only with ``refinement="boxes"``: with ``"tagging"`` it is
       silently ignored.
     - None: no refinement, one level. Boxes of one level must not overlap, and
       the number of levels is the number of keys plus one.
   * - ``max_nbr_levels``
     - Total number of levels, level 0 included. You must give it with
       ``"tagging"``. With ``"boxes"`` it is deduced from ``refinement_boxes``
       and giving it is an error.
     - 1 (no refinement). With ``"tagging"`` and no value, a warning is
       printed. PHARE refuses more than 10 levels or fewer than 1 when the
       simulator starts.
   * - ``tagging_threshold``
     - ``"tagging"`` only. A cell is tagged when a measure of how fast the
       magnetic field changes around it is above this number. In 1D: how much
       the 5-cell average of By and Bz changes from one cell to the next,
       relative to 1 + that average. In 2D and 3D: the jump of Bx, By or Bz
       over two cells, relative to 1 + the jump over one cell. Values are in
       code units, so the threshold is relative to a field of order 1. Lower
       value, more cells tagged, larger fine regions.
     - 0.1. No check in Python.
   * - ``tag_buffer``
     - ``"tagging"`` only. Number of cells added around each tagged cell
       before the fine patches are built. Larger value, wider fine regions
       around the tagged cells.
     - 1. No check in Python.
   * - ``clustering``
     - How the tagged cells are grouped into patches: ``"tile"`` (SAMRAI's
       tile clustering) or ``"berger"`` (the Berger-Rigoutsos algorithm).
     - ``"tile"``. Anything else raises ``ValueError``.
   * - ``nesting_buffer``
     - Each fine level must stay at least this many cells of the level below
       away from the edge of that level (a number, or one per direction). In
       ``"boxes"`` mode PHARE checks your boxes against it. For level 0 it
       means: keep refined boxes away from the box edge.
     - 1. Must be at least 1 and at most half of ``smallest_patch_size``.
       Giving it without refinement only prints a warning and has no effect.
   * - ``smallest_patch_size``
     - Smallest allowed patch, in cells (a number or one per direction). A
       box in ``refinement_boxes``, once refined, must be at least this large.
     - 3. Must be larger than 2 and not larger than ``cells``. With
       ``interp_order=1`` (the default), 4 and 5 are refused.
   * - ``largest_patch_size``
     - Largest patch, in cells, in each direction (a number or one per
       direction).
     - None (no limit). Must be at least 1, at least ``smallest_patch_size``,
       not larger than ``cells``.
   * - ``refined_particle_nbr``
     - How many fine particles replace one coarse particle where a fine
       region is created. Only some values are allowed.
     - First allowed value for your dimension and ``interp_order``
       (1D, default ``interp_order=1``: 2, allowed 2 and 3; order 2:
       allowed 2, 3, 4).
   * - ``refinement_order``
     - Order of the operator that fills the fine fields from the coarse ones.
     - 2, the only value allowed.

You cannot set the refinement ratio: giving ``refinement_ratio`` is an error.
For how the patches are shared between processes, see :doc:`load_balancing`.
To see what was tagged, write the ``"tags"`` quantity with the diagnostics
(:doc:`diagnostics`: ``ph.MetaDiagnostics(quantity="tags", write_timestamps=...)``).

For the theory, see :doc:`/numerics/amr`.


Common mistakes
---------------

* ``ValueError: Box(Box([0],[9])) is incompatible with coarser boxes([Box([0],[39])]) and nest_buffer([1])``
  — your box touches the edge of the level below, or sticks out of it. Move it
  at least ``nesting_buffer`` cells inside. On a 40-cell box with the default,
  level-0 boxes must lie within cells 1 to 38.
* The same message for an ``"L1"`` box: ``"L1"`` boxes use level-1 numbers.
  With ``"L0": (10,),(19,)``, level 1 has cells 20 to 39, and an ``"L1"``
  box must lie inside cells 21 to 38. Do not give level-0 numbers there.
* ``ValueError: Invalid box incompatible with smallest_patch_size`` — the box,
  once refined, is smaller than ``smallest_patch_size`` (3 by default). A box
  of 1 level-0 cell gives 2 fine cells and is refused.
* ``ValueError: Error: Box(...) overlaps with Box(...)`` — two boxes of the
  same level share cells. Merge them or separate them.
* ``ValueError: Error: invalid arguments - max_nbr_levels`` — you gave
  ``max_nbr_levels`` with ``refinement="boxes"``. Remove it, the number of
  levels comes from ``refinement_boxes``. (Keep it with ``"tagging"``.)
* ``ValueError: Error: invalid arguments - refinement_ratio`` — the ratio is
  fixed at 2. Remove it.
* ``WARNING, 'max_nbr_levels' is not set, defaulting to 1`` — with
  ``refinement="tagging"`` you forgot ``max_nbr_levels``. The run goes on with
  no refinement at all.
* ``refinement="tagging"`` with a ``refinement_boxes`` — no message, the boxes
  are ignored.
* A typo such as ``refinement="tag"`` — no message from Python. With
  ``max_nbr_levels`` above 1 the run crashes at start (segmentation fault);
  with ``max_nbr_levels=1`` it runs without refinement. Write ``"tagging"``
  exactly.
* ``ValueError: Error: nesting_buffer([0]) must be at least 1 on a simulation with refinement levels``
  — use 1 or more.
* ``ValueError: Error: clustering type is not supported, supported types are ['berger', 'tile']``
  — use one of these two.
* ``WARNING, nesting_buffer(2) has no effect on a simulation with no refinement levels, ignoring``
  — you gave ``nesting_buffer`` but no refinement.

.. seealso:: How this works inside: :doc:`/dev/one_time_step`.
