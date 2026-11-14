===================
The patch hierarchy
===================

**Goal:** understand the object you get when you read PHARE output with
``pharesee``: its levels, its patches, the ghost cells, and how to loop over
all of it and get the numbers you want.

.. contents:: On this page
   :local:
   :depth: 1


Minimal example
---------------

PHARE refines the mesh where needed: a coarse level covers the whole box, and
finer levels cover only some regions. Each level is cut into *patches*,
rectangular pieces of grid. The data you read back keeps this structure,
and it is called a *patch hierarchy* (class ``PatchHierarchy``).

First, a small 1D run with two levels: level 0 covers the box, and level 1
covers a part of it (the cells 10 to 19 of level 0, given with
``refinement_boxes``). Save it as ``amr_run.py`` and run it with ``python3
amr_run.py`` (a few seconds). It is a small setup in the style of :doc:`/getting_started/first_simulation`,
with a refined region added; see :doc:`/user_guide/amr` for the
refinement options.

.. code-block:: python

   import numpy as np
   import pyphare.pharein as ph
   from pyphare.simulator.simulator import Simulator

   sim = ph.Simulation(
       cells=40, dl=0.5, time_step=0.01, final_time=0.1,
       refinement_boxes={"L0": {"B0": [(10,), (19,)]}},
       diag_options={"format": "phareh5",
                     "options": {"dir": "amr_run", "mode": "overwrite"}},
   )
   def density(x): return 1.0
   def bx(x): return 1.0
   def by(x): return 0.1*np.cos(2*np.pi*x/20)
   def bz(x): return 0.0
   def vx(x): return 0.0
   def vth(x): return 0.1
   ph.MaxwellianFluidModel(
       bx=bx, by=by, bz=bz,
       protons={"charge": 1, "density": density,
                "vbulkx": vx, "vbulky": vx, "vbulkz": vx,
                "vthx": vth, "vthy": vth, "vthz": vth,
                "nbr_part_per_cell": 100},
   )
   ph.ElectronModel(closure="isothermal", Te=0.01)
   ph.ElectromagDiagnostics(quantity="B", write_timestamps=[0.0, 0.1])
   Simulator(sim).run()

Then read the output. :doc:`reading_data` explains ``Run`` and ``GetB``; here
we look at what comes back.

.. code-block:: python

   from pyphare.pharesee.run import Run

   run = Run("amr_run")                          # folder written by the simulation
   print(run.times("B"))                         # output times of B: [0.  0.1]

   B = run.GetB(0.0, all_primal=False)           # a PatchHierarchy, at t = 0
   print(B.quantities())                         # ['Bx', 'By', 'Bz']
   print(B.levelNbrs(), B.finest_level())        # [0, 1] 1

   for ilvl, level in B.levels(0.0).items():
       print("level", ilvl, "cell width", level.cell_width)
       for patch in level.patches:
           pd = patch["By"]
           g = pd.ghosts_nbr[0]
           print("  patch", patch.id, "box", patch.box, "origin", patch.origin,
                 "dl", patch.dl)
           print("    By:", pd.centerings[0], pd.dataset.shape[0], "points,",
                 g, "ghost cells per side")
           x = pd.x                              # positions of the points of pd.dataset
           x_domain = x[g:-g]                    # the same without the ghosts
           by_domain = pd.dataset[g:-g]
           print("    x from", x[0], "to", x[-1], "- without ghosts", x_domain[0], x_domain[-1])
           print("    pd[box] has", pd[patch.box].shape[0], "points")

   fine = B.level(1).patches[0]
   print("level 1 cells", fine.box.lower, "to", fine.box.upper,
         "= level 0 cells", fine.box.lower // 2, "to", fine.box.upper // 2)

   merged = run.GetB(0.0, merged=True)
   by, (x,) = merged["By"]
   print(by(7.3))

On this run it prints:

.. code-block:: text

   [0.  0.1]
   ['Bx', 'By', 'Bz']
   [0, 1] 1
   level 0 cell width [0.5]
     patch p0#0 box Box([0],[39]) origin [0.] dl [0.5]
       By: dual 44 points, 2 ghost cells per side
       x from -0.75 to 20.75 - without ghosts 0.25 19.75
       pd[box] has 40 points
   level 1 cell width [0.25]
     patch p0#0 box Box([20],[39]) origin [5.] dl [0.25]
       By: dual 24 points, 2 ghost cells per side
       x from 4.625 to 10.375 - without ghosts 5.125 9.875
       pd[box] has 20 points
   level 1 cells [20] to [39] = level 0 cells [10] to [19]
   -0.06791865

(The patch name ``p0#0`` and the very last number may differ on your machine.
With several processes, each level has several patches.)

``all_primal=False`` asks for the fields as PHARE stores them (see "Primal and
dual" below). Without it, ``GetB`` returns every component on primal points
and under the names ``"x"``, ``"y"``, ``"z"`` (see "Common mistakes").


What the object contains
------------------------

The picture shows a hierarchy like a larger version of the one above: two
levels, two patches on each. Each patch is drawn as a row of cells.

.. plot::
   :caption: A two-level 1D hierarchy (schematic, not read from a run). Level 0
             has cells of width 1, level 1 has cells of width 1/2. Each patch has
             a box: its first and last cell, numbered in the index space of its
             own level. The hatched cells are ghost cells (2 per side here); the
             rest are the patch's own cells.

   import matplotlib.pyplot as plt
   from matplotlib.patches import Rectangle

   ng = 2  # ghost cells per side (interpolation order 1)
   # (level, lower, upper) in that level's own cell indices
   patches = [(0, 0, 19), (0, 20, 39), (1, 20, 39), (1, 56, 71)]
   dl = {0: 1.0, 1: 0.5}
   rowy = {0: [3.0, 2.0], 1: [0.9, -0.1]}
   colors = {0: "#9ecae1", 1: "#fdae6b"}

   fig, ax = plt.subplots(figsize=(5.2, 4.2))
   count = {0: 0, 1: 0}
   for lvl, lo, up in patches:
       y = rowy[lvl][count[lvl]]
       count[lvl] += 1
       w = dl[lvl]
       for i in range(lo - ng, up + ng + 1):
           ghost = i < lo or i > up
           ax.add_patch(Rectangle((i * w, y), w, 0.7,
                                  facecolor="0.85" if ghost else colors[lvl],
                                  edgecolor="0.35" if lvl == 0 else "0.45", lw=0.5,
                                  hatch="///" if ghost else None))
       ax.text((lo + up + 1) / 2 * w, y + 0.8, f"box = [{lo}, {up}]", ha="center", va="bottom", fontsize=13)
       ax.text(-3.2, y + 0.35, f"level {lvl}", ha="right", va="center", fontsize=13)

   ax.axvline(0, color="k", lw=1)
   ax.axvline(40, color="k", lw=1)
   ax.set_xlabel("x (physical position)", fontsize=13)
   ax.text(43.5, -0.6, "hatched = ghost cells", ha="right", va="center", fontsize=13)
   ax.set_xlim(-3.5, 43.5)
   ax.set_ylim(-0.9, 4.6)
   ax.set_yticks([])
   ax.set_xticks([0, 10, 20, 30, 40])
   ax.tick_params(labelsize=13)
   for s in ("left", "right", "top"):
       ax.spines[s].set_visible(False)
   ax.text(0.0, 4.55, "domain edge", ha="center", va="bottom", fontsize=13, clip_on=False)
   ax.text(40.0, 4.55, "domain edge", ha="center", va="bottom", fontsize=13, clip_on=False)
   fig.tight_layout()


**Levels.** Level 0 covers the whole box. Level 1 has cells twice as small,
level 2 twice as small again, and so on: the cell width on level ``n`` is the
level 0 width divided by ``2**n``. In the run above, level 1 covers only part
of the box.

**Patches.** A level is a list of patches. A patch is a rectangle of cells
(a segment in 1D).

**The box.** ``patch.box`` gives the first and last cell of the patch (both
included) as integer indices. These indices belong to the index space of the
patch's level: on every level, cell 0 is the first cell at the lower edge of
the simulation box, and the whole box of level ``n`` has ``2**n`` times more
cells than level 0. ``hierarchy.level_domain_box(n)`` gives the whole box of
level ``n``. In the run above, level 0 has the box ``[0, 39]`` (40
cells) and level 1 has the box ``[0, 79]``.

**Same cell on two levels.** One cell of level ``n`` is made of two cells of
level ``n + 1`` (in 1D; four in 2D). Cell ``i`` of level 0 is therefore the
cells ``2 i`` and ``2 i + 1`` of level 1. This is why the level 1 patch above
has the box ``[20, 39]`` while we asked for the level 0 cells 10 to 19. The
same factor applies from level 1 to level 2. **Never compare indices from two
different levels without converting them.**

**Origin and cell width.** ``patch.origin`` is the position of the lower edge
of the patch's first cell, and ``patch.dl`` is the cell width on that level.
Check: level 1 patch above, ``5.0 = 20 * 0.25``.

**Time.** A hierarchy can hold several output times. The time is not stored in
the patch: it is the key you use to get the levels, ``hierarchy.levels(t)``.
Without a time, you get the first one.

**Quantities (patch data).** Each patch carries one *patch data* per
quantity (here ``Bx``, ``By``, ``Bz``). Fields are stored in a ``FieldData``
object: ``pd.dataset`` is the array of numbers (for data read from a file it
is read on demand; ``pd.dataset[:]`` gives a NumPy array), ``pd.x`` the
positions of its points.

**Ghost cells.** A patch data holds the numbers of the patch's own cells
**and** a few ghost cells on each side. PHARE fills the ghost cells during the run, because the
numerical scheme needs values around the cells it updates at the edge of a patch. Both ``pd.dataset`` and ``pd.x`` include
them; ``pd.ghosts_nbr`` is how many per side and per direction (2 here;
this number comes from the file). To get only the patch's own points, drop
``g = pd.ghosts_nbr[0]`` points at each end: ``pd.dataset[g:-g]``. How the
ghost cells are filled is described in :doc:`/dev/ghost_cells`.

**Primal and dual.** PHARE does not store all quantities at the same
places in a cell. A quantity lives either on the *primal* points (the cell
edges, on the grid lines) or on the *dual* points (the cell centers), in each
direction. In 1D, ``Bx`` is primal and ``By`` and ``Bz`` are dual.
``pd.centerings`` gives ``"primal"`` or ``"dual"`` for each direction. This changes the number of points. For a box of ``N`` cells and ``g``
ghost cells per side:

* dual: ``N + 2 g`` points. In the run above, level 0 ``By``: 40 + 4 = 44.
* primal: ``N + 2 g + 1`` points (the last edge closes the last cell). Level
  0 ``Bx``: 40 + 4 + 1 = 45.

and the dual points are shifted by half a cell: for ``By`` on level 0 the
first point is at ``-0.75 = 0 - 2 * 0.5 + 0.25``.


How to
------

**Loop over levels and patches.**

.. code-block:: python

   for ilvl, level in B.levels(0.0).items():     # {level number: PatchLevel}
       for patch in level.patches:                # or: for patch in level
           ...

**Choose a time.** ``B = run.GetB(0.0)`` reads that one time. To read several
times into one hierarchy, give a list, ``run.GetB([0.0, 0.1], all_primal=False)``,
then ``B.times()`` lists them and ``B.levels(0.1)`` selects one. Note
that ``B.times()`` returns the times as text (``'0.0000000000'``), while
``run.times("B")`` returns numbers. Both ``0.1`` and the text can be
used to select.

**x coordinates of a patch field.** ``pd.x`` (and ``pd.y``, ``pd.z`` in 2D
and 3D) with the ghosts, ``pd.x[g:-g]`` without. They have as many entries as
``pd.dataset`` has along that direction.

**The own points of a patch, as a box.** ``pd[patch.box]`` returns the
patch's own points as a NumPy array (no ghosts): ``N`` points for
dual, ``N + 1`` for primal. You can pass any ``Box``, in the indices of the
patch's level: you get the points inside it, ghost points included if the box
reaches into them, and an empty array if the box does not touch the patch data.

**The value at a position, from the finest level that has it.** Use
``merged=True``, as in the minimal example: ``run.GetB(t, merged=True)`` gives,
for each component, a function ``f`` and the grid coordinates, and ``f(x)``
takes the finest data available at ``x`` (nearest grid point by default,
``interp="nearest"``; see
:doc:`reading_data`).

**The finest data as flat arrays.** ``flat_finest_field(B, "By", time=0.0)``
from ``pyphare.pharesee.hierarchy.hierarchy_utils`` returns two 1D arrays
(values, positions) made from the finest patches, where coarse points under a finer
level are dropped. It keeps one ghost point at each patch end and the points
are not sorted by position; it is what ``merged=True`` uses.
In 2D the values and positions are flattened, in 3D it is not implemented.


Main attributes and methods
---------------------------

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Name
     - What it gives
   * - ``hierarchy.levels(time=None)``
     - dict ``{level number: PatchLevel}`` at that time (first time if none)
   * - ``hierarchy.level(n, time=None)``
     - the ``PatchLevel`` of level ``n``
   * - ``hierarchy.levelNbrs(time=None)``, ``levelNbr(time=None)``
     - list of level numbers, and their count
   * - ``hierarchy.finest_level(time=None)``
     - the highest level number
   * - ``hierarchy.times()``
     - the output times, as an array of text
   * - ``hierarchy.quantities()``
     - names of the patch data, e.g. ``['Bx', 'By', 'Bz']``
   * - ``hierarchy.domain_box``, ``hierarchy.ndim``
     - box of level 0 (cells), and the dimension
   * - ``hierarchy.level_domain_box(n)``
     - the whole box of level ``n``, in the index space of level ``n``
   * - ``level.patches``, ``level[i]``, ``for patch in level``
     - the patches of the level (a list), one by index, or all in turn
   * - ``level.cell_width``
     - cell width of the level (from its first patch)
   * - ``level.level_number``
     - the number of the level
   * - ``level.level_range()``
     - lowest and highest ``x`` of the first quantity's points, ghosts included (1D)
   * - ``patch.box``
     - first and last own cell of the patch, in its level's indices
   * - ``patch.origin``, ``patch.dl``
     - position of the lower edge of the first cell; cell width
   * - ``patch.id``
     - the patch name (a text, here ``p0#0``)
   * - ``patch["By"]``
     - the patch data named ``By``
   * - ``patch.patch_datas``
     - dict ``{name: patch data}``
   * - ``pd.dataset``
     - the numbers, ghosts included
   * - ``pd.x``, ``pd.y``, ``pd.z``
     - positions of the points of ``dataset`` along each direction
   * - ``pd.centerings``
     - ``"primal"`` or ``"dual"``, per direction
   * - ``pd.ghosts_nbr``
     - ghost cells per side, per direction
   * - ``pd.ghost_box``
     - the box grown by the ghosts (cells)
   * - ``pd[box]``
     - the points inside ``box``; ``pd[patch.box]`` gives the own points
   * - ``pd.box``, ``pd.origin``, ``pd.dl``
     - same as the patch's


Common mistakes
---------------

* ``KeyError: No patchdata for key: Bx`` on a hierarchy from ``GetB(t)`` — the
  default ``GetB`` renames the components to ``x``, ``y``, ``z`` and puts them
  on primal points. Use ``patch["x"]``, or ``GetB(t, all_primal=False)`` to keep
  ``Bx``, ``By``, ``Bz``.
* Ghost points are ``NaN`` with the default ``GetB`` — on that default, the
  ghost points are not computed and hold ``NaN`` (4 points in the run above,
  2 per side). Take ``pd.dataset[g:-g]`` (with ``g = pd.ghosts_nbr[0]``) for the own points. With
  ``all_primal=False`` the ghost points hold numbers.
* An array has more points than the box has cells — it includes the ghost
  points, and one more for primal data. Use ``pd.dataset[g:-g]`` or
  ``pd[patch.box]``, and read ``pd.centerings``.
* Indexing the array with the cell index of the box — ``pd.dataset[0]`` is the
  first ghost point, not the first cell of the box. The first own point is
  ``pd.dataset[g]``.
* Comparing or mixing indices of two levels — the box of a level 1 patch is in
  level 1 indices: multiply level 0 indices by 2 (or divide by 2) before
  comparing, as in the example above.
* Summing or averaging ``pd.dataset`` over all patches of all levels —
  the arrays include ghost points, and a refined region has data on both the
  coarse and the fine level, so it is counted twice. Use the data without
  ghosts and, where levels overlap, only the finest one (``merged=True`` does
  this).

.. seealso:: How ghost cells are filled: :doc:`/dev/ghost_cells`. To plot a
   hierarchy, see :doc:`plotting_fields`.
