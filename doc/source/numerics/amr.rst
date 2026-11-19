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

:term:`Refinement <refinement>` computes fine values from coarse ones, at one given time.
PHARE's rule is **second order** (the only order available: ``refinement_order`` must be
2) and the refinement ratio is 2. It works one direction at a time:

* where a fine node sits on a coarse node, it takes the coarse value;
* where it sits halfway between two coarse nodes, it takes their average;
* along a direction where values sit at cell centres, each coarse cell
  :math:`I` gives two fine values, :math:`c_I \mp (c_{I+1}-c_{I-1})/8`: a
  linear profile through the cell.

In 2D and 3D the rule is the product of the 1D rules.

.. plot::
   :include-source: False
   :alt: One dimension. Cell-centred values: three coarse values give two fine values. Node-centred values: a fine node copies a coarse node or averages two.

   import matplotlib.pyplot as plt
   from phare_draw import C

   c = [1.0, 1.2, 1.6]                         # coarse cell values I-1, I, I+1
   s = (c[2] - c[0]) / 8
   left, right = c[1] - s, c[1] + s
   assert abs(left - 1.125) < 1e-12 and abs(right - 1.275) < 1e-12 and abs((left + right) / 2 - c[1]) < 1e-12
   n = [0.1, 0.3]                              # two coarse nodes
   mid = (n[0] + n[1]) / 2
   assert abs(mid - 0.2) < 1e-12
   fig, axes = plt.subplots(2, 1, figsize=(3.4, 5.0))
   ax = axes[0]
   for k, v in enumerate(c):
       ax.add_patch(plt.Rectangle((k, 1.0), 1, 0.45, fc="white", ec=C["L0"], lw=1.8))
       ax.text(k + 0.5, 1.22, f"{v}", ha="center", va="center", fontsize=9)
   for k, v in enumerate((left, right)):
       ax.add_patch(plt.Rectangle((1 + 0.5 * k, 0.2), 0.5, 0.45, fc=C["patch"], ec=C["L1"], lw=1.8))
       ax.text(1.25 + 0.5 * k, 0.42, f"{v:.3f}", ha="center", va="center", fontsize=9, color=C["L1"])
   ax.text(1.5, 1.62, "cell-centred: two children of the middle cell", ha="center", fontsize=8)
   ax.text(1.5, -0.05, "their average is 1.2, the coarse value", ha="center", fontsize=8, color=C["grey"])
   ax.set_xlim(-0.1, 3.1); ax.set_ylim(-0.2, 1.8); ax.set_axis_off()
   ax = axes[1]
   for x, v in zip((0, 1), n):
       ax.plot([x], [1.0], "s", mfc="white", mec=C["L0"], mew=1.8, ms=10)
       ax.text(x, 1.22, f"{v}", ha="center", fontsize=9)
   for x, v, t in [(0, n[0], "copy"), (0.5, mid, "average"), (1, n[1], "copy")]:
       ax.plot([x], [0.2], "o", color=C["L1"], ms=8)
       ax.text(x, 0.0, f"{v:.1f}", ha="center", va="top", fontsize=9, color=C["L1"])
       ax.text(x, 0.55, t, ha="center", fontsize=8, color=C["grey"])
   ax.text(0.5, 1.62, "node-centred: copy or average", ha="center", fontsize=8)
   ax.set_xlim(-0.3, 1.3); ax.set_ylim(-0.4, 1.8); ax.set_axis_off()
   fig.tight_layout()

What the rule keeps:

* the two children of a cell average back to the coarse value, also in 2D and 3D;
* a field that varies linearly in space is reproduced exactly (the children lie on the
  line through the coarse value with the slope given by the neighbours), which is what
  second order means here.

The same rules serve hybrid and MHD. MHD refines its conservative variables (density,
momentum, total energy), which sit at cell centres, with the cell-centre rule: a refined
coarse cell keeps its mass, momentum and energy (for the cells filled from the coarse level).

The magnetic field is refined in two stages so that its divergence is kept: each component
is first refined with the rule above, then the fine faces inside each coarse cell are
corrected so that every fine cell has the divergence of the coarse cell it comes from. The
two stages, a drawing with numbers and what this guarantees across levels are on
:doc:`div_b`. The correction follows Balsara, Samantaray and Subramanian, *Efficient
WENO-based prolongation strategies for divergence-preserving vector fields*, Commun. Appl.
Math. Comput., `doi:10.1007/s42967-021-00182-x <https://doi.org/10.1007/s42967-021-00182-x>`_.
How the code does it is on :doc:`/dev/field_refinement`.


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
* the ion charge density, which lives on the grid nodes: each coarse node takes the value
  of the fine node on top of it;
* the ion current on the electric-field positions, which Ohm's law reads
  (:doc:`hybrid_ohm`): like the electric field;
* the magnetic field is not averaged. The coarse level redoes its Faraday
  update over the coarse step with the fine electric field, averaged in time
  over the fine sub-steps, where the fine level exists (this is called
  :term:`refluxing <reflux>`). The magnetic flux through a coarse face then
  changes by the sum of the changes through the fine faces it contains.

**MHD** does not coarsen its state. At the same moment, the fine fluxes and electric field
summed over the sub-steps replace the coarse fluxes on the coarse faces they cover,
and the coarse level redoes its step with them (again refluxing). Details:
:doc:`/dev/refluxing` and :doc:`/dev/coarsening`.

What refluxing keeps, with a fixed time step:

* **MHD:** the mass, momentum and energy that leave a coarse cell through a face of the
  border are exactly what the fine cells on the other side receive, because the coarse
  face is given the time average of the fine fluxes (a quarter of each of the four
  sub-step fluxes) and, in 2D and 3D, the mean over the 2 (2D) or 4 (3D) fine faces it
  covers.
* **Hybrid and MHD:** the magnetic flux through a coarse face of the border changes by the
  sum of the changes through the fine faces it contains, so the divergence of B is kept
  across the border (:doc:`div_b`).

.. plot::
   :include-source: False
   :alt: Flux through a coarse-fine face over one coarse step: four fine sub-step fluxes and the single coarse flux, before and after refluxing.
   :caption: Illustrative fluxes through a face of the coarse–fine border over one coarse step. After refluxing, the area under the grey line (coarse flux times the coarse step) equals the sum of the areas of the four blue bars (fine fluxes times a quarter of the step).

   import numpy as np
   import matplotlib.pyplot as plt
   from phare_draw import C

   F_fine = np.array([1.10, 1.20, 1.30, 1.40])    # illustrative fluxes of the four fine sub-steps
   F_coarse = 1.00                                # illustrative flux of the coarse step
   mean = F_fine.mean()
   assert abs(mean - 1.25) < 1e-12
   assert abs(np.sum(F_fine * 0.25) - mean * 1.0) < 1e-12    # time integrals over the coarse step agree
   fig, axes = plt.subplots(2, 1, figsize=(3.4, 4.6), sharex=True, sharey=True)
   for ax, title in zip(axes, ["before: the coarse level used its own flux", "after: it uses the mean of the fine fluxes"]):
       for k, f in enumerate(F_fine):
           ax.add_patch(plt.Rectangle((0.25 * k, 0), 0.25, f, fc=C["patch"], ec=C["L1"], lw=1.5))
           ax.text(0.25 * k + 0.125, f / 2, f"{f:.2f}", ha="center", va="center", fontsize=8, color=C["L1"])
       ax.set_title(title, fontsize=9, loc="left")
       ax.set_xlim(0, 1); ax.set_ylim(0, 1.7)
       ax.set_yticks([]); ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
       ax.set_xticklabels(["t", "", "", "", "t + Δt"], fontsize=8)
       for s in ("top", "right", "left"):
           ax.spines[s].set_visible(False)
   axes[0].plot([0, 1], [F_coarse] * 2, color=C["L0"], lw=3)
   axes[0].text(0.5, 1.58, f"coarse flux {F_coarse:.2f} (grey)", ha="center", fontsize=8, color=C["L0"])
   axes[1].plot([0, 1], [mean] * 2, color=C["L0"], lw=3)
   axes[1].text(0.5, 1.58, f"coarse flux {mean:.2f} = mean (grey)", ha="center", fontsize=8, color=C["L0"])
   fig.text(0.5, 0.01, "illustrative numbers", ha="center", fontsize=7, color=C["grey"])
   fig.tight_layout(rect=(0, 0.03, 1, 1))

With an adaptive time step in MHD the number of fine sub-steps can differ from four and the
mean is no longer the time average: see the known issue on :doc:`/dev/refluxing`.


Fields at level boundaries
--------------------------

Each patch stores a few layers of :term:`ghost cells <ghost cell>` around
its own cells, so that derivatives can be computed at its edge. Where the ghost cells
overlap another patch of the same level, their values are copied from it
(:term:`patch ghosts <patch ghost>`). Where they lie outside the level, there is no data
of the same level (:term:`level ghosts <level ghost>`). They are filled in two steps, in
this order:

#. the coarse values are **interpolated in time**, linearly, between the values saved at
   the start of the coarse step and the values at its end:
   :math:`(1-\alpha)\,\text{old}+\alpha\,\text{new}`;
#. the result is **refined in space** with the rules above.

The time interpolation is exact for a field that varies linearly in time over the coarse
step, and its error is second order in the time step.

.. plot::
   :include-source: false
   :alt: One coarse step, four fine sub-steps, the two saved coarse times, and the times at which B and the averaged E are filled in sub-step 2 (hybrid).
   :caption: Hybrid, sub-step 2 (shaded): the coarse level holds an *old* value (saved at the start of its step) and a *new* one (already computed for the end of the step); the averaged electric field (E) is filled at the middle of the sub-step and B at its end.

   from phare_draw import *
   fig, (a,) = strip(1, aspect=False, height=3.0, phone=True)
   amr_clock(a, side_labels=True, shade=[2],
             coarse_dots=[(0, "old"), (1, "new")],
             marks=[(0.375, "E"), (0.5, "B")])
   finish(fig, [a], align=False)

* **Hybrid:** each field is filled at its own time inside the sub-step, so each has its own
  weight :math:`\alpha`. B, Bpred and E are filled at the end of the sub-step:
  :math:`\alpha=0.25, 0.5, 0.75, 1` for the four sub-steps. The averaged electric field is
  filled at the middle of the sub-step: :math:`\alpha=0.125, 0.375, 0.625, 0.875`.
* **MHD:** the ghosts of density, momentum, total energy and B are filled at the time of
  each Runge–Kutta stage inside each sub-step, so the weights depend on the integrator
  (:doc:`mhd_time_integration`).

How the code does it: :doc:`/dev/time_interpolation`. The ghost B, and the cells next to the
level border where it is not divergence-free, are discussed on :doc:`div_b`.

The order of the whole scheme with refinement is measured by the 2D AMR convergence tests of
the repository (``WENOZ``, MHD, an Alfvén wave and a whistler wave): the measured order must
lie between 1.70 and 2.25 (:doc:`/validation/mhd`).


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

Ohm's law reads the ion charge density and ion current one node beyond the level
(:doc:`hybrid_ohm`), where the level ghost particles do not reach in full. There, both are
filled like the fields (hybrid): refined in space from the coarser level and interpolated in
time between its values at the start and at the end of the coarse step. Where the
level ghost particles are complete, the fine level keeps its own particle values.


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
