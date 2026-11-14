.. _dev-gridlayout:

=======================================
GridLayout: centerings and index spaces
=======================================

.. container:: breadcrumb

   You are here: :doc:`/dev/index` › GridLayout

**In one sentence:** a ``GridLayout`` describes one :term:`patch`'s grid: where each
quantity sits in a cell (:term:`primal` or :term:`dual` per direction, the
:term:`Yee layout`), how many :term:`ghost <ghost cell>` values surround the patch, and
how an array index maps to an AMR index and to a position. User view: :doc:`/numerics/grid`.

The picture
-----------

.. plot::
   :include-source: False
   :caption: 1 · One direction, a patch of 4 cells starting at AMR cell 8, 2 ghosts
      (interpolation order 1). Numbers above the markers: local array index; bottom row: AMR
      cell index. 2 · ``project<Bz, Ex>`` in 2D: the Ex value at the circled edge is the
      average of the two Bz around it.

   from phare_draw import *

   fig, ax = strip(2, titles=["Local and AMR indices, primal and dual",
                              "project<Bz, Ex>(Bz, i, j) in 2D"])
   a = ax[0]
   lo, hi, g = 8, 12, 2
   for x0 in (lo - g, hi):
       a.add_patch(__import__("matplotlib").patches.Rectangle((x0, -1.5), g, 3.0, fc=C["ghost"], ec="none", zorder=0))
   a.plot([lo - g, hi + g], [0, 0], color=C["grey"], lw=1)
   for k, x in enumerate(range(lo - g, hi + g + 1)):
       ghost = x < lo or x > hi
       mark(a, "primal", x, 0.6, color=C["grey"] if ghost else C["E"])
       text(a, x, 1.05, str(k), fs=12)
   for k in range(hi - lo + 2 * g):
       x = lo - g + k + 0.5
       ghost = x < lo or x > hi
       mark(a, "dual", x, -0.5, color=C["grey"] if ghost else C["B"])
       text(a, x, -0.95, str(k), fs=12)
   for c in range(lo - g, hi + g):
       text(a, c + 0.5, -1.35, str(c), fs=12, color=C["grey"])
   text(a, lo - g - 0.2, 0.6, "primal", ha="right", color=C["E"])
   text(a, lo - g - 0.2, -0.5, "dual", ha="right", color=C["B"])
   text(a, lo - g - 0.2, -1.35, "AMR cell", ha="right", color=C["grey"])
   text(a, (lo + hi) / 2, 1.45, "domain", bold=True)
   text(a, lo - g / 2, 1.45, "ghosts", color=C["ghost_edge"])
   text(a, hi + g / 2, 1.45, "ghosts", color=C["ghost_edge"])
   a.update_datalim([(lo - g - 2.6, 0)])

   a = ax[1]
   grid(a, (0, 0, 2, 2))
   for y in (0.5, 1.5):
       mark(a, "dual", 0.5, y, color=C["B"])
       arrow(a, (0.5, y), (0.5, 1 + (y - 1) * 0.35), color=C["grey"], lw=1.4, ms=10)
       text(a, 0.82, y, "Bz", color=C["B"])
   a.add_patch(__import__("matplotlib").patches.Circle((0.5, 1), 0.13, fc=C["hi"], ec=C["E"], lw=2, zorder=4))
   text(a, 0.12, 1.0, "Ex", color=C["E"], bold=True, ha="right")
   finish(fig, ax)

* **One layout per patch**, built on demand by ``layoutFromPatch`` from the SAMRAI patch
  box, mesh size and origin (dl = 1, origin = 0 for SAMRAI's temporary patches that have no
  geometry) (:src:`src/amr/resources_manager/amr_utils.hpp#L169-L209`).
* **Centering table per model**, a compile-time array indexed by the quantity's enum
  (``initLayoutCentering_``, :src:`src/core/data/grid/impl/yee/gridlayout_hybrid_yee.hpp#L52-L98`).
  MHD differs: ρ, ρv, P, Etot at cell centres
  (:src:`src/core/data/grid/impl/yee/gridlayout_mhd_yee.hpp#L60-L80`).
* **Index tables** per centering and direction: first/last domain index, last ghost index;
  the first ghost index is 0 (:src:`src/core/data/grid/gridlayout.hpp#L1131-L1238`).
* **Loop boxes**: ``evalOnBox`` (domain), ``evalOnGhostBox``, ``evalOnShrinkedGhostBox``,
  ``evalOnBiggerBox`` (:src:`src/core/data/grid/gridlayout.hpp#L973-L997`).
* **Two index spaces**: local (array position, 0 = first ghost) and AMR (cell index on the
  level, shared by all its patches), ``localToAMR`` / ``AMRToLocal``
  (:src:`src/core/data/grid/gridlayout.hpp#L582-L660`). Patch exchanges speak AMR,
  numerical loops speak local.

The tricks
----------

Same number of ghosts for primal and dual
   **Problem.** Primal has one more value than dual; different ghost counts would need
   offsets in every loop.
   **Fix.** Both get ``field_ghost_width`` ghosts, so they start at the same local index; only
   the end differs, ``start + N`` vs ``start + N − 1``
   (``nDNbrGhosts``, :src:`src/core/data/grid/gridlayout.hpp#L533-L537`; ``dualOffset_``,
   :src:`src/core/data/grid/gridlayout.hpp#L1087`).

Dual index i is half a cell after primal index i
   **Problem.** Derivatives and averages need a fixed neighbour convention.
   **Fix.** Position :math:`i\,\Delta x` (primal), :math:`(i+\tfrac12)\Delta x` (dual), from
   the AMR index (``fieldNodeCoordinates``, :src:`src/core/data/grid/gridlayout.hpp#L450-L476`).
   So primal → dual uses :math:`i, i+1`; dual → primal uses :math:`i-1, i`
   (``Partial``, :src:`src/core/numerics/differential_operators.hpp#L47-L82`).

Averaging stencils derived from centerings
   **Problem.** 58 hand-written averaging tables, one per pair of quantities.
   **Fix.** ``GridLayout::project<Source, Destination>`` builds the stencil at compile time:
   two-point average in each direction where the centerings differ, product across
   directions (default order 2; the machinery also has orders 4 and 6)
   (``recenteringStencil``, :src:`src/core/data/grid/gridlayoutdefs.hpp#L259-L275`).
   Debug builds check the field's centering
   (:src:`src/core/data/grid/gridlayout.hpp#L662-L672`). :pr:`1499` (stencils derived from
   centerings), :pr:`1495` (centering lookup table).

Ghost width from the numerics
   **Problem.** Ghosts must cover the widest stencil and whole coarse cells (ratio 2).
   **Fix.** Hybrid: 2, 4, 4 for interpolation order 1, 2, 3
   (:src:`src/core/models/options/hybrid_options.hpp#L24-L25`). MHD: reconstruction width
   + 1, rounded up to even, i.e. 2 (``Constant``) or 4
   (:src:`src/core/utilities/ghost_width_calculator.hpp#L47-L68`).

Reading a field at another centering, outside solvers
   **Problem.** Diagnostics and tagging want every field at one place, without a copy.
   **Fix.** ``RecenteredField``: a read-only view that averages when read; it compares
   source and destination centerings on every read, so it is not for solver loops
   (:src:`src/core/data/field/recentered_field.hpp#L1-L98`). :pr:`1486` (read at another
   centering).

Go deeper
---------

* :doc:`/dev/fields` — the arrays the layout describes.
* :doc:`/dev/ghost_cells` — how ghost values are filled.

Related
-------

* :doc:`/dev/electromagnetic` · :doc:`/dev/mhd_solver` · :doc:`/numerics/grid`
