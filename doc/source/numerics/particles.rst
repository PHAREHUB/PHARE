========================================
Particles: shapes, pusher, interpolation
========================================

**Goal:** how PHARE's macro-particles represent the ions, move, and talk to
the grid, and which options control this.

Minimal example
---------------

The particle options of a hybrid run (commented lines), on top of
:doc:`/getting_started/first_simulation`:

.. code-block:: python

   import pyphare.pharein as ph
   from pyphare.simulator.simulator import Simulator

   sim = ph.Simulation(
       cells=200, dl=0.5, time_step=0.01, final_time=1.0,
       interp_order=2,                   # particle shape: 1, 2 or 3 (default 1)
       particle_pusher="modified_boris", # the only pusher (default)
       diag_options={"format": "phareh5",
                     "options": {"dir": "particles_run", "mode": "overwrite"}},
   )

   ph.MaxwellianFluidModel(
       bx=1.0,
       protons={"charge": 1, "density": 1.0, "vthx": 0.1, "vthy": 0.1, "vthz": 0.1,
                "nbr_part_per_cell": 200,   # macro-particles per cell (default 100)
                "init": {"seed": 1337}},    # same particles at every run
   )
   ph.ElectronModel(closure="isothermal", Te=0.01)

   Simulator(sim).run()

What a macro-particle is
------------------------

* A :term:`macro-particle` stands for many ions of one :term:`population`:
  position, velocity, *weight*; charge and mass are the population's.
* At start, ``nbr_part_per_cell`` particles per cell where the density is above
  ``density_cut_off``: Maxwellian velocities, uniform positions in the cell.
* Weight = local density / ``nbr_part_per_cell``: a density per particle.

Shapes: particles are small clouds
----------------------------------

For the grid a particle is a cloud centred on its position. ``interp_order``
(:term:`interpolation order`) picks the B-spline shape: order 1 (2 cells wide),
2 (3 cells) or 3 (4 cells). Higher orders: smoother moments, less noise, more
cost, wider ghost regions.

.. plot::
   :include-source: False
   :caption: Particle at x = 0; distances in cells. Each shape has unit area.

   import numpy as np
   from phare_draw import *
   x = np.linspace(-2.5, 2.5, 1001)
   u = np.abs(x)
   shapes = {
       "order 1": np.clip(1 - u, 0, None),
       "order 2": np.where(u < 0.5, 0.75 - u**2, np.where(u < 1.5, 0.5 * (1.5 - u)**2, 0.0)),
       "order 3": np.where(u < 1, 2/3 - u**2 + u**3 / 2, np.where(u < 2, (2 - u)**3 / 6, 0.0)),
   }
   fig, (a,) = strip(1, ["The three particle shapes"], height=2.4, aspect=False)
   for (label, sh), c, ty in zip(shapes.items(), (C["L1"], C["L2"], C["ion"]), (0.95, 0.78, 0.62)):
       a.plot(x, 2.5 * sh, lw=2.4, color=c)
       text(a, 1.75, 2.5 * ty, label, color=c, bold=True, ha="left")
   for v in range(-2, 3):
       mark(a, "primal", v, 0, color=C["grey"], size=6); text(a, v, -0.3, str(v), color=C["grey"])
   finish(fig, [a])

The same shape is used both ways (the usual particle-in-cell choice):

* **particles → grid** (*deposition*): density += weight × shape, flux +=
  weight × velocity × shape, on the nodes the cloud covers.
* **grid → particles** (*interpolation*): E and B at the particle are the
  shape-weighted average of nearby values, each component at its place on the
  staggered grid (:doc:`/numerics/grid`).

The pusher: modified Boris
--------------------------

Particles move with a Boris scheme (the :term:`pusher`), the standard leap-frog
for charged particles. One time step Δt:

.. plot::
   :include-source: False
   :caption: Top: the position moves in two half drifts, with the fields read
      at the middle. Bottom: the velocity update in the plane normal to B (B out
      of the page): half kick by E, rotation around B, half kick by E.

   import numpy as np
   from phare_draw import *
   fig, (a1, a2) = strip(2, ["Position: two half drifts", "Velocity: kick, rotate, kick"], height=2.6)
   grid(a1, (0, 0, 6, 2), color=C["light"])
   p0, p1, p2 = (0.7, 0.6), (2.6, 1.0), (4.9, 1.5)
   particles(a1, [p0], color=C["grey"], size=50); particles(a1, [p1], color=C["L2"], size=50)
   particles(a1, [p2], color=C["ion"], size=50)
   arrow(a1, p0, p1, color=C["grey"], lw=1.6, ms=10, text="vⁿ Δt/2")
   arrow(a1, p1, p2, color=C["ion"], lw=1.6, ms=10, text="vⁿ⁺¹ Δt/2")
   text(a1, 0.7, 0.25, "n"); text(a1, 2.6, 0.62, "n+½: read E, B", color=C["L2"]); text(a1, 4.9, 1.15, "n+1")
   o = np.array([0.0, 0.0]); kick = np.array([0.45, 0.0])
   vn = np.array([1.2, 0.5]); vm = vn + kick
   th = -np.deg2rad(45); R = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
   vp = R @ vm; v1 = vp + kick
   arrow(a2, o, vn, color=C["grey"], lw=1.8, ms=10); text(a2, *(vn * 0.5 + [-0.1, 0.22]), "vⁿ")
   arrow(a2, vn, vm, color=C["E"], lw=1.8, ms=10, text="(q/m) E Δt/2")
   ang = np.linspace(np.arctan2(vm[1], vm[0]), np.arctan2(vp[1], vp[0]), 30)
   r = np.linalg.norm(vm)
   a2.plot(r * np.cos(ang), r * np.sin(ang), color=C["B"], lw=1.6, ls="--")
   text(a2, 1.95, 0.05, "rotation\naround B", color=C["B"], ha="left")
   arrow(a2, vp, v1, color=C["E"], lw=1.8, ms=10, text="(q/m) E Δt/2", text_off=(0, -0.25))
   mark(a2, "node", *vm, color=C["B"]); text(a2, *(vm + [0.12, -0.05]), "v⁻", color=C["B"], ha="left")
   mark(a2, "node", *vp, color=C["B"]); text(a2, *(vp + [-0.2, 0.05]), "v⁺", color=C["B"])
   arrow(a2, o, v1, color=C["ion"], lw=2.2, ms=10); text(a2, *(v1 + [0.2, -0.05]), "vⁿ⁺¹", ha="left")
   mark(a2, "node", 0, 0, color=C["ink"])
   a2.plot([-0.3, 2.6], [-0.5, 0.9], alpha=0)
   finish(fig, [a1, a2])

The fields are the average of those at the start of the step and the predicted
ones at its end (:doc:`/numerics/time_integration`). Particles are pushed twice
per step: once to predict the ion moments, once for real.

A particle may move about **2 cells per half step** at most (a bit less forward,
more backward); beyond that the run stops (see below). Safe rule:
:math:`v_\text{max}\,\Delta t / \Delta x < 2` (under one cell per half step).

Refinement: splitting particles
-------------------------------

In a refined region each coarse particle becomes ``refined_particle_nbr`` finer
ones around it, same velocity, nearly the same density. Possible values (bold:
default; compiled in a standard build):

.. list-table::
   :header-rows: 1

   * - Dimension
     - Order 1
     - Order 2
     - Order 3
   * - 1D
     - **2**, 3
     - **2**, 3, 4
     - **2**, 3, 4, 5
   * - 2D
     - **4**, 5, 8, 9
     - **4**, 5, 8, 9, 16
     - **4**, 5, 8, 9, 25
   * - 3D
     - **6**, 12, 27
     - **6**, 12, 27
     - **6**, 12, 27

Options
-------

.. list-table::
   :header-rows: 1
   :widths: 25 20 40 15

   * - Option
     - Where
     - Meaning
     - Default
   * - ``interp_order``
     - ``Simulation``
     - particle shape order: 1, 2 or 3
     - 1
   * - ``particle_pusher``
     - ``Simulation``
     - particle pusher; only ``"modified_boris"`` exists
     - ``"modified_boris"``
   * - ``refined_particle_nbr``
     - ``Simulation``
     - particles made from one coarse particle when refining (table above)
     - 2 / 4 / 6 in 1D / 2D / 3D
   * - ``nbr_part_per_cell``
     - population dict
     - particles loaded per cell
     - 100
   * - ``density_cut_off``
     - population dict
     - cells with a lower density get no particles
     - 1e-16
   * - ``init``
     - population dict
     - ``{"seed": <int>}``: reproducible particle loading, whatever the patch
       layout or number of MPI processes
     - not seeded

Common mistakes
---------------

* **The run aborts with only "MPI_ABORT was invoked" in the terminal.** Look in
  the log of each :term:`MPI` process, ``.log/0.out`` (``.log/<rank>.out``), in
  the folder you launched the run from (with the default ``PHARE_LOG``, see
  :doc:`/user_guide/environment_variables`). If it contains ``PHARE ERROR: ... PrePush Particle moved 2 cells``
  (or ``PostPush``), a particle went too far in half a step: reduce
  ``time_step``, or check the velocities of your initial condition.
* **Choosing a ``refined_particle_nbr`` that is not in the table.** In 3D, the
  input check accepts more values than the table, but a standard build has no
  compiled simulator for them, so the run cannot start.
* **Noisy moments.** Increase ``nbr_part_per_cell`` or ``interp_order``; both
  cost run time.

.. seealso:: How this works inside: :doc:`/dev/particles_ions`.
