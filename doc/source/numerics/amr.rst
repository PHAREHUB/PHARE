========================
Adaptive Mesh Refinement
========================

:term:`Adaptive mesh refinement <adaptive mesh refinement>` (AMR) makes the
mesh finer only where the solution needs it. This page explains the algorithm
PHARE uses: how the refined mesh is organized, how it is advanced in time,
and how data moves between coarse and fine meshes. How to set it up is on
:doc:`/user_guide/amr`.

.. contents:: On this page
   :local:
   :depth: 1


Patch based approach
--------------------

The mesh is a hierarchy of :term:`levels <level>`. Level 0 covers the whole
domain. Each finer level covers part of the level below, with cells twice
smaller in every direction: the :term:`refinement ratio` is always 2. Up to
10 levels are allowed.

A level is a set of rectangular :term:`patches <patch>`. Each patch is a
small grid with its own arrays, and patches are shared between MPI processes.
A fine level must lie inside the level below it, at least ``nesting_buffer``
coarse cells away from its edge (*proper nesting*), so that every fine cell
near the edge of the fine level still has coarse cells around it to take
boundary values from.

.. plot::
   :include-source: false
   :alt: A 2D hierarchy: level 0 made of four patches, level 1 made of two patches inside it.

   from phare_draw import *
   from matplotlib.patches import Rectangle

   fig, (ax,) = strip(1, ["a two-level hierarchy"], phone=True)
   for (x, y) in [(0, 0), (16, 0), (0, 8), (16, 8)]:
       ax.add_patch(Rectangle((x, y), 16, 8, fc="none", lw=2, ec=C["L0"]))
   for (x, y, w, h) in [(6, 4, 10, 6), (16, 4, 8, 8)]:
       ax.add_patch(Rectangle((x, y), w, h, fc=C["patch"], ec=C["L1"], lw=2, alpha=0.9))
   text(ax, 8, 12, "level 0", color=C["L0"])
   text(ax, 15, 7.2, "level 1:\n2 patches", color=C["L1"])
   ax.plot([-1, 33], [-1, 17], alpha=0)
   finish(fig, [ax])

*Schematic of a two-level hierarchy. Level 1 (blue) lies inside level 0
(grey), away from the domain edge.*


Recursive time integration
--------------------------

The hybrid equations carry whistler waves, whose stable time step scales like
:math:`\Delta x^2`. A level twice finer therefore uses a time step four times
smaller, and takes 4 sub-steps for each step of the level below, as long as the level-0 step does not change, which is always the case in hybrid. With an adaptive step (MHD only) the number can differ (:doc:`/dev/subcycling`).

Levels are advanced coarsest first, recursively: a level makes one step, then
every finer level makes its sub-steps to reach the same time, then the levels
are synchronized.

.. plot::
   :include-source: false
   :alt: Order in which three levels are advanced during one coarse time step.

   from phare_draw import *
   from matplotlib.patches import Rectangle

   order = []
   def recurse(level, t0, t1):
       order.append((level, t0, t1))
       if level < 2:
           dt = (t1 - t0) / 4
           for k in range(4):
               recurse(level + 1, t0 + k * dt, t0 + (k + 1) * dt)
   recurse(0, 0.0, 16.0)

   fig, (ax,) = strip(1, ["one coarse step, three levels"])
   col = [C["L0"], C["L1"], C["L2"]]
   for i, (lvl, t0, t1) in enumerate(order, start=1):
       y = 2 - lvl
       ax.add_patch(Rectangle((t0, y * 1.5), t1 - t0, 1.1, fc=col[lvl], ec="white", lw=1.5))
       text(ax, (t0 + t1) / 2, y * 1.5 + 0.55, str(i), color="white", bold=True,
            fs=13 if lvl < 2 else 11, halo=False)
   for lvl in range(3):
       text(ax, -0.4, (2 - lvl) * 1.5 + 0.55, f"level {lvl}", ha="right")
   for k, lab in enumerate([r"$t$", r"$t+\Delta t/4$", r"$t+\Delta t/2$", r"$t+3\Delta t/4$", r"$t+\Delta t$"]):
       text(ax, 4 * k, -0.45, lab)
   ax.plot([-3.2, 16.5], [-0.8, 4.2], alpha=0)
   finish(fig, [ax])

*One coarse step with three levels: numbers give the order of the 21 level
advances.*

Because the coarse level goes first, a fine level always has the coarse
solution at the beginning and at the end of its sub-steps. The values it
needs at its border during a sub-step are interpolated in time between these
two coarse states (see "Fields at level boundaries" below). How the code
implements this is described in :doc:`/dev/subcycling`.


Regridding
----------

Regridding is the process of changing the geometry of an existing mesh level
in the AMR hierarchy, by changing the number and geometry of some patches of
that level. This typically occurs when the needs for refinement evolve with
the solution over time.

Concretely, regridding means removing an entire level from the AMR hierarchy
and replacing it by a new one. The so-called "new" level is initialized by
copying values from the old level wherever they overlap, and by refining
values from coarser levels where they do not.

Technically, regridding involves different methods depending on whether
particle data or field data is refined. The specific operation performed also
depends on the type of field being refined: for instance, the magnetic field
refinement must preserve the divergence-free character of the field.

In PHARE, a level's finer levels are regridded every 2 steps of that level.
With fixed refinement boxes the new levels have the same boxes; with tagging
they follow the tagged cells (see :ref:`tagging-for-refinement`).


Field refinement
----------------

:term:`Refinement <refinement>` computes fine values from coarse ones. PHARE
uses a second-order rule in each direction:

* where a fine node sits on a coarse node, it takes the coarse value;
* where it sits halfway between two coarse nodes, it takes their average;
* along a direction where values sit at cell centres, each coarse cell
  :math:`I` gives two fine values, :math:`c_I \mp (c_{I+1}-c_{I-1})/8`: a
  linear profile through the cell.

In 2D and 3D the rule is the product of the 1D rules.

The magnetic field is refined in two stages so that its divergence stays
zero: each component is first refined with the rule above, then the fine
faces inside each coarse cell are corrected so that every fine cell has the
divergence of the coarse cell it comes from, which is zero. The correction
follows Balsara, Samantaray and Subramanian, *Efficient WENO-based
prolongation strategies for divergence-preserving vector fields*, Commun.
Appl. Math. Comput., `doi:10.1007/s42967-021-00182-x
<https://doi.org/10.1007/s42967-021-00182-x>`_.

The same rules serve hybrid and MHD. MHD refines its conservative variables
(density, momentum, total energy), which sit at cell centres, with the cell-centre
rule (the third one): its two children average back to the coarse value, so a
refined coarse cell keeps its mass, momentum and energy. How the code does it is on
:doc:`/dev/field_refinement`.


Particle refinement
-------------------

Where a new fine level is created, each coarse macro-particle is replaced by
several lighter fine particles placed around it (:term:`particle splitting`).
The number of fine particles per coarse particle is
``refined_particle_nbr`` (:term:`refined particle number`). The splitting
patterns are chosen so that the fine particles deposit nearly the same
density as the coarse particle they replace. Details are on
:doc:`/numerics/particles`.


Field coarsening
----------------

Coarsening, also sometimes called "restriction", is the process of projecting
some quantities of a given level onto the region it overlaps on the next
coarser level (which ones depends on the model, see below). Coarsening is used so that the solution on the coarse level
"feels" the fine solution in the overlapped region, which is assumed to be of
better quality.

When a fine level has caught up in time with its coarser level, **hybrid**
coarsens:

* the electric field: each coarse edge gets the average of the 2 fine edges
  along it or, for components that sit on nodes, the value of the fine node on
  top of it;
* the ion charge density and charge flux, which live on the grid nodes: each
  coarse node takes the value of the fine node on top of it;
* the magnetic field is not averaged. The coarse level redoes its Faraday
  update over the coarse step with the fine electric field, averaged in time
  over the fine sub-steps, where the fine level exists (this is called
  :term:`refluxing <reflux>`). The magnetic flux through a coarse face then
  changes by the sum of the changes through the fine faces it contains.

**MHD** does not coarsen its state. At the same moment, the fine fluxes and electric field
summed over the sub-steps replace the coarse fluxes on the coarse faces they cover,
and the coarse level redoes its step with them (again refluxing). Details:
:doc:`/dev/refluxing` and :doc:`/dev/coarsening`.


Fields at level boundaries
--------------------------

Each patch stores a few layers of :term:`ghost cells <ghost cell>` around
its own cells, so that derivatives can be computed at its edge. Where the
ghost cells overlap another patch of the same level, their values are copied
from it (:term:`patch ghosts <patch ghost>`). Where they lie outside the
level, there is no data of the same level: values are refined from the
coarser level and interpolated linearly in time between the two coarse
states that bracket the fine sub-step (:term:`level ghosts <level ghost>`).
Each field is filled at its own time within the sub-step, so its time weight
differs; see :doc:`/dev/time_interpolation`.


Particles at level boundaries
-----------------------------

Plasma flows across the border of a refined level, which is not a physical
boundary. To let particles enter the fine level and to complete the moments
near its border, coarse particles near the border are split into a thin layer
of fine particles just outside the level, the
:term:`level ghost particles <level ghost particle>`. At each coarse step two
sets exist, one from the coarse particles at the start of the coarse step and
one from the coarse particles at its end; their contributions to the moments
are weighted by the time of the fine sub-step.


.. _tagging-for-refinement:

Tagging for refinement
-----------------------

With ``refinement="tagging"``, PHARE marks (:term:`tags <tagging>`) the cells
where the magnetic field varies quickly, before each regrid:

* in 1D, how much the 5-cell average of :math:`B_y` and :math:`B_z` changes
  from one cell to the next, relative to 1 + that average;
* in 2D and 3D, the jump of :math:`B_x`, :math:`B_y` or :math:`B_z` over two
  cells, relative to 1 + the jump over one cell.

A cell is tagged when this measure exceeds ``tagging_threshold``. The tagged
region is then grown by ``tag_buffer`` cells in every direction, so that a
moving structure stays inside the refined region until the next regrid.


.. _clustering:

Clustering
----------

Tagged cells must be grouped into rectangular patches. PHARE offers two
algorithms, chosen with ``clustering``.

.. _tiling:

Tiling
~~~~~~

``"tile"`` (the default) uses SAMRAI's tile clustering.

.. _berger-rigoutsos:

Berger-Rigoutsos
~~~~~~~~~~~~~~~~

``"berger"`` uses the Berger-Rigoutsos algorithm, which covers the tagged
cells with boxes by recursively splitting a bounding box where few cells are
tagged. PHARE limits these boxes to 10 cells in each direction.
