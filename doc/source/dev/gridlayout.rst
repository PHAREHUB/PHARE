.. _dev-gridlayout:

==========================================
GridLayout: centerings and index spaces
==========================================

.. container:: breadcrumb

   You are here: :doc:`/dev/index` › GridLayout

**In one sentence:** a ``GridLayout`` describes one patch's grid: where each physical
quantity sits in a cell (its *centering*), how many ghost nodes surround the patch,
and how an array index turns into an AMR index and a position in space.

Why it exists
-------------

PHARE stores each field component in its own array, on a staggered grid (the *Yee
layout*): B components sit on cell faces, E and J on cell edges, and the hybrid ion
moments on cell corners (*nodes*). Every numerical operator needs to know, for every
array it reads, where its values are. Writing that by hand in every loop is how index
bugs are born. ``GridLayout`` answers these questions in one place, and most answers are
computed at compile time.

The user-side description of the grid (pictures of the cell, what lives where) is on
:doc:`/numerics/grid`.

The picture
-----------

One direction, one patch of 4 cells, interpolation order 1 (so 2 ghost nodes on each
side). Each centering has its own array:

.. plot::
   :include-source: False
   :caption: Primal values (circles) sit on cell boundaries (nodes), dual values (squares) at cell
      centres. Both arrays start with the same number of ghost nodes (grey band), so
      their first domain value has the same *local* index, 2. The primal array has one
      more value than the dual one. Below the axis: the *AMR* cell index of the same
      places on the level, here with the patch starting at AMR cell 8.

   import matplotlib.pyplot as plt

   fig, ax = plt.subplots(figsize=(6.5, 3.4))
   lo, hi = 8, 12            # domain cells 8..11 (AMR index)
   g = 2                     # ghost width
   ax.axvspan(lo - g, lo, color="0.85", lw=0)
   ax.axvspan(hi, hi + g, color="0.85", lw=0)
   for x in range(lo - g, hi + g + 1):
       ax.plot([x, x], [-0.15, 0.15], color="0.5", lw=1)
   ax.plot([lo - g, hi + g], [0, 0], color="0.3", lw=1)
   # primal: local 0..8 at x = 6..14
   for k, x in enumerate(range(lo - g, hi + g + 1)):
       ghost = x < lo or x > hi
       ax.plot(x, 0.55, "o", ms=10, color="0.6" if ghost else "#2a7ab9")
       ax.text(x, 0.85, str(k), ha="center", fontsize=11)
   # dual: local 0..7 at x = 6.5..13.5
   for k in range(hi - lo + 2 * g):
       x = lo - g + k + 0.5
       ghost = x < lo or x > hi
       ax.plot(x, -0.55, "s", ms=9, color="0.6" if ghost else "#e07a1f")
       ax.text(x, -0.95, str(k), ha="center", fontsize=11)
   for c in range(lo - g, hi + g):
       ax.text(c + 0.5, -1.45, str(c), ha="center", fontsize=11, color="0.25")
   ax.text(lo - g - 0.3, 0.85, "primal, local", ha="right", fontsize=11)
   ax.text(lo - g - 0.3, -0.95, "dual, local", ha="right", fontsize=11)
   ax.text(lo - g - 0.3, -1.45, "AMR cell", ha="right", fontsize=11, color="0.25")
   ax.text((lo + hi) / 2, 1.25, "domain (4 cells)", ha="center", fontsize=12)
   ax.text(lo - g / 2, 1.25, "ghosts", ha="center", fontsize=11)
   ax.text(hi + g / 2, 1.25, "ghosts", ha="center", fontsize=11)
   ax.set_xlim(2.2, hi + g + 0.4)
   ax.set_ylim(-1.7, 1.5)
   ax.axis("off")
   fig.tight_layout()

How it works
------------

#. **One layout per patch, built on demand.** ``layoutFromPatch`` reads the SAMRAI
   patch box (its AMR cells), the mesh size and the origin, and builds a ``GridLayout``
   (:src:`src/amr/resources_manager/amr_utils.hpp#L169-L200`). It is a small value object;
   solvers build one per patch inside their loops.
#. **A centering table per model.** For each physical quantity and direction, the table
   says primal or dual. The hybrid table (``initLayoutCentering_``,
   :src:`src/core/data/grid/impl/yee/gridlayout_hybrid_yee.hpp#L52-L98`) puts B on faces,
   E and J on edges, and every ion moment on nodes. The MHD table
   (:src:`src/core/data/grid/impl/yee/gridlayout_mhd_yee.hpp`) differs: density,
   momentum, pressure and total energy are at cell centres (dual in every direction).
#. **Index tables.** From the ghost width and the number of cells, the layout fills, per
   centering and direction, the first and last *physical* (domain) index and the last
   *ghost* index (``initPhysicalStart_``, ``initPhysicalEnd_``, ``initGhostEnd_``,
   :src:`src/core/data/grid/gridlayout.hpp#L1131-L1238`). The first ghost index is always 0.
#. **Boxes to loop on.** Operators do not write loops by hand. They ask for a box and a
   function: ``evalOnBox`` (domain values), ``evalOnGhostBox`` (all values),
   ``evalOnShrinkedGhostBox`` and ``evalOnBiggerBox`` (in between)
   (:src:`src/core/data/grid/gridlayout.hpp#L973-L997`). The box depends on the field's
   centering, so a primal field loops over one more value than a dual one.
#. **Two index spaces.** The *local* index is the position in the patch's array; it starts
   at 0 on the first ghost. The *AMR* index is the cell index on the whole level, shared by
   all patches of that level. ``localToAMR`` and ``AMRToLocal`` convert between them
   (:src:`src/core/data/grid/gridlayout.hpp#L582-L660`). Data exchange between patches
   (ghost filling, copies) speaks AMR indices; numerical loops speak local indices.

The tricks
----------

Primal and dual get the same number of ghosts
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** A primal array has one more value than a dual one over the same cells. If
the ghost counts also differed, every loop would need per-centering offsets.

**Fix.** Both centerings get ``field_ghost_width`` ghosts (``nDNbrGhosts``,
:src:`src/core/data/grid/gridlayout.hpp#L533-L537`; ``initPhysicalStart_``,
:src:`src/core/data/grid/gridlayout.hpp#L1131-L1153`). So the first domain value has the
same local index for both, and only the *end* differs: with :math:`N` cells, primal
domain values run from ``start`` to ``start + N``, dual ones from ``start`` to
``start + N - 1`` (``dualOffset_``, :src:`src/core/data/grid/gridlayout.hpp#L1087`).

Dual index :math:`i` is half a cell *after* primal index :math:`i`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** A derivative or an average between a primal and a dual array needs a fixed
convention for which neighbours to take.

**Fix.** The value of index :math:`i` is at :math:`i\,\Delta x` if primal and at
:math:`(i + \tfrac12)\,\Delta x` if dual (``fieldNodeCoordinates``,
:src:`src/core/data/grid/gridlayout.hpp#L450-L476`, computed from the AMR index). Hence:

* a derivative of a primal field at dual index :math:`i` uses primal :math:`i` and
  :math:`i+1`; a derivative of a dual field at primal index :math:`i` uses dual
  :math:`i-1` and :math:`i` (``Partial``,
  :src:`src/core/numerics/differential_operators.hpp#L47-L82`);
* averaging to the other centering uses the same pairs.

Projections are derived from the centerings
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** Ohm's law needs, for example, :math:`B_z` on the :math:`E_x` edges, where
:math:`B_z` does not live. Each such pair used to have its own hand-written averaging
table.

**Fix.** ``GridLayout::project<Source, Destination>(field, index)`` builds the averaging
stencil at compile time from the two centerings (``recenteringStencil``,
:src:`src/core/data/grid/gridlayoutdefs.hpp#L259-L275`). In every direction where they
differ, it averages the two neighbours; across several directions it takes the product.
In 2D, :math:`B_z` is dual–dual and :math:`E_x` is dual–primal, so

.. math::

   B_z \text{ on } E_x(i, j) = \tfrac12\big(B_z(i, j-1) + B_z(i, j)\big).

Debug builds check that the field passed really has the source centering
(:src:`src/core/data/grid/gridlayout.hpp#L662-L672`).

**History.** #1499 (stencils derived from centerings).

Centerings are a compile-time table
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** The centering of a quantity is asked inside hot loops.

**Fix.** It is a ``constexpr`` array indexed by the quantity's enum value
(``centering``, :src:`src/core/data/grid/impl/yee/gridlayout_hybrid_yee.hpp#L127-L144`):
no switch, no runtime cost, and a missing or duplicate entry fails to compile.

**History.** #1495 (centering lookup table).

The ghost width comes from the numerics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** Ghosts must be wide enough for the widest stencil that reads them, and a
ghost layer must cover whole coarse cells so that refinement from the coarser level
(ratio 2) fills it cleanly.

**Fix.** The ghost width is a compile-time option of the model:

* hybrid: 2, 4, 4 for interpolation orders 1, 2, 3
  (:src:`src/core/models/options/hybrid_options.hpp#L24-L25`), set by the particle
  shape;
* MHD: the reconstruction stencil width plus one, rounded up to an even number
  (``nbrGhostsFromReconstruction``,
  :src:`src/core/utilities/ghost_width_calculator.hpp#L47-L68`): 2 for ``Constant``, 4 for
  the others. The extra one is there because J is computed on the ghost box minus one
  layer and the reconstruction still reads it.

Reading a field at another centering, outside the solver
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** Diagnostics and tagging want every field at one place (for example cell
centres), without allocating a copy.

**Fix.** ``RecenteredField`` is a read-only view that averages on the fly
(:src:`src/core/data/field/recentered_field.hpp`). It checks centerings at run time, so it
is meant for whole-field reads, not for solver loops, which use ``project``.

**History.** #1486 (read at another centering).

Go deeper
---------

* :doc:`/dev/fields`: the arrays the layout describes.
* :doc:`/dev/ghost_cells`: how ghost values are filled.

Related
-------

* :doc:`/dev/electromagnetic`: the main users of ``project`` and the derivatives.
* :doc:`/dev/mhd_solver`: the MHD layout in use.
* :doc:`/numerics/grid`: the same grid from the user's side.
