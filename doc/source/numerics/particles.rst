========================================
Particles: shapes, pusher, interpolation
========================================

**Goal:** know how PHARE's macro-particles represent the ions, how they move,
how they exchange information with the grid, and which input options control
this.

Minimal example
---------------

The particle options of a hybrid run, on top of the
:doc:`/getting_started/first_simulation` script (only the lines that matter here
are commented):

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

A :term:`macro-particle` stands for many ions of one :term:`population`. It has a position, a
velocity and a *weight*. The charge and mass are those of its population.

At the start, PHARE loads ``nbr_part_per_cell`` particles in every cell where
the density is above ``density_cut_off``. Velocities are drawn from a
Maxwellian with the bulk and thermal velocities you give, positions uniformly in
the cell. All particles of a cell get the same weight, the density there divided
by ``nbr_part_per_cell``: the weight is a density per particle.

Shapes: particles are small clouds
----------------------------------

A particle is not a point for the grid: it is a cloud of a given *shape*,
centred on its position. The ``interp_order`` option (:term:`interpolation order`) chooses the shape: the
B-spline of order 1 (a triangle, 2 cells wide), 2 (3 cells) or 3 (4 cells).
Higher orders give smoother moments and less noise, at a higher cost per
particle, and need wider ghost regions around each patch.

.. plot::
   :include-source: False
   :caption: The three particle shapes, for a particle at x = 0, in cells.
      Each has unit area.

   import numpy as np
   import matplotlib.pyplot as plt

   x = np.linspace(-2.5, 2.5, 1001)
   u = np.abs(x)
   shapes = {
       "order 1": np.clip(1 - u, 0, None),
       "order 2": np.where(u < 0.5, 0.75 - u**2, np.where(u < 1.5, 0.5 * (1.5 - u)**2, 0.0)),
       "order 3": np.where(u < 1, 2/3 - u**2 + u**3 / 2, np.where(u < 2, (2 - u)**3 / 6, 0.0)),
   }
   fig, ax = plt.subplots(figsize=(6, 3.2))
   for (label, s), c in zip(shapes.items(), ("#2a7ab9", "#f4a259", "#3a9a5b")):
       ax.plot(x, s, lw=2.2, label=label, color=c)
   ax.set_xlabel("distance from the particle (cells)", fontsize=12)
   ax.tick_params(labelsize=11)
   ax.legend(fontsize=11, frameon=False)
   ax.grid(axis="x", color="0.9")
   for side in ("top", "right"):
       ax.spines[side].set_visible(False)
   fig.tight_layout()

The same shape is used both ways:

* **particles → grid** (*deposition*): each particle adds weight × shape to the
  density, and weight × velocity × shape to the flux, on the grid nodes its cloud
  covers. Summing over particles gives the ion moments.
* **grid → particles** (*interpolation*): the electric and magnetic fields at
  a particle are the shape-weighted average of the field values on the nodes
  around it, taking into account where each field component sits on the
  staggered grid (:doc:`/numerics/grid`).

Using the same shape both ways is the usual choice in particle-in-cell codes:
it keeps deposition and interpolation consistent and limits spurious
self-forces.

The pusher: modified Boris
--------------------------

Particles move with a Boris scheme (the :term:`pusher`), the standard leap-frog for charged
particles. Over one time step :math:`\Delta t`:

.. mermaid::

   flowchart TD
     A["position at n"] -->|"half drift: x += vⁿ Δt/2"| B["position at n+1/2"]
     B --> C["interpolate E, B<br/>at the particle"]
     C --> D["velocity n → n+1:<br/>half kick by E, rotation by B,<br/>half kick by E"]
     D -->|"half drift: x += vⁿ⁺¹ Δt/2"| E["position at n+1"]

The fields used are the average of the fields at the start of the step and the
predicted fields at its end (see :doc:`/numerics/time_integration`). In each
time step, particles are pushed twice: once to predict the ion moments at the
new time, once for real.

A particle may move about **2 cells in half a time step** at most (a bit less
forward, a bit more backward, depending on where it starts in its cell). Beyond
that, the run stops with an error (see below). This bounds the time step by the fastest
particles: with :math:`v_\text{max}\,\Delta t / \Delta x < 2`, no particle moves
more than one cell in half a step, which is always safe.

Refinement: splitting particles
-------------------------------

When a region is refined, coarse particles are replaced by
``refined_particle_nbr`` finer particles placed around them, with the same
velocity, so that they deposit nearly the same density on the finer grid. The
default number depends on the dimension; other values are possible:

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

Bold: default. These are the combinations compiled in a standard build.

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
