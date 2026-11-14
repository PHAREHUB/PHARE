:orphan:

Drawing kit (for authors)
=========================

Not part of the reader-facing guide. Shows every helper of ``doc/source/_ext/phare_draw.py`` (matplotlib, drawn at build time)
and the ``.. flow::`` directive. Copy a block, change the numbers.

Yee cell, patch with ghosts, 1D values and particles
----------------------------------------------------

.. plot::
   :include-source: True
   :caption: B on the faces (magenta) and E on the edges (blue) of one cell; a patch with its ghost cells (dashed tint); primal (dots on cell edges) and dual (squares at cell centres) values with particles inside the cells.

   from phare_draw import *
   fig, (a, b, c) = strip(3, ["One Yee cell: B left, E right", "A patch (6 x 3 cells) with 2 ghost cells", "Values and particles along x"], height=2.9)
   yee_pair(a)
   grid(b, (-2, -2, 8, 5)); patch(b, (0, 0, 6, 3), ghosts=2, label="patch 0")
   ruler(c, 0, 6)
   particles(c, [(0.3, .3), (0.8, .35), (2.2, .3), (2.5, .35), (4.6, .3)], color=C["ion"], size=22)
   legend(c, [("primal (node)", C["rho"], "primal"), ("dual (centre)", C["B"], "dual"), ("particle", C["ion"], "particle")], loc="upper center", bbox=(0.5, 1.0), ncol=3)
   finish(fig, [a, b, c])

Level and boxes
---------------

.. plot::
   :include-source: True
   :caption: A refined region (ratio 2, fine cells drawn) over a coarse grid; boxes with arrows for layers.

   from phare_draw import *
   fig, (a, b) = strip(2, ["Coarse grid, one refined region", "Layers"], height=3.0)
   grid(a, (0, 0, 8, 4)); level(a, (2, 1, 5, 3), "level 1", ratio=2)
   arrow(a, (7.3, 2.0), (5.15, 2.0), color=C["L1"]); text(a, 6.4, 2.45, "fine cells, dx/2", color=C["L1"])
   boxes(b, ["Python\ninput, run control", "C++ engine\nsimulator, amr", "core\nfields, particles"], arrows=[(0, 1, "calls"), (1, 2, "uses")])
   finish(fig, [a, b])

Flow (clickable boxes)
----------------------

.. flow::

   Start | the box is a link | /dev/index
   then
   Second box | small text

Phone-sized figure
------------------

Few elements per panel? Use ``phone=True``: the figure is 3.4 in wide (340 px), so it is seen whole on a phone.
Dense figures keep the default 6.5 in width and scroll sideways on a phone.

.. plot::
   :include-source: True
   :caption: A phone-sized strip: three steps, each panel simple enough to read at 340 px.

   from phare_draw import *
   fig, (a, b, c) = strip(3, ["Patch", "With its ghost cells", "Refined x2"], height=2.0, phone=True)
   grid(a, (0, 0, 4, 3)); patch(a, (0, 0, 4, 3), label="patch")
   grid(b, (-1, -1, 5, 4)); patch(b, (0, 0, 4, 3), ghosts=1)
   grid(c, (0, 0, 4, 3)); level(c, (1, 1, 3, 2), ratio=2)
   finish(fig, [a, b, c])
