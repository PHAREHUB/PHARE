==========
MHD inputs
==========

**Goal:** set up and run a magnetohydrodynamics (:term:`MHD`) simulation: choose
the numerical scheme, describe the initial fluid, and read the result back.

The example is a one-dimensional shock tube (two uniform states separated by a
jump). It runs in a few seconds.

.. contents:: On this page
   :local:
   :depth: 1


Minimal example
---------------

Work in an empty folder and make sure Python finds ``pyphare`` and your build
folder (see :doc:`/getting_started/first_simulation`). Save this as
``shock.py``:

.. code-block:: python

   import numpy as np
   import pyphare.pharein as ph
   from pyphare.simulator.simulator import Simulator
   from pyphare.pharesee.run import Run

   ph.NO_GUI()

   L, ncells = 400.0, 400
   x0 = L / 2                      # position of the jump
   final_time = 40.0

   # 1. The simulation: an MHD run, and the numerical scheme to use
   sim = ph.Simulation(
       time_step=0.2,
       final_time=final_time,
       cells=ncells,
       dl=L / ncells,
       model_options=["MHDModel"],     # this is an MHD run
       gamma=5.0 / 3.0,                # adiabatic index
       reconstruction="Linear",
       limiter="VanLeer",
       riemann="Rusanov",
       mhd_timestepper="TVDRK2",
       diag_options={"format": "phareh5",
                     "options": {"dir": "mhd_run", "mode": "overwrite"}},
   )

   # 2. The initial fluid: left and right states
   ph.MHDModel(
       density=lambda x: np.where(x < x0, 1.0, 0.125),
       vx=0.0,
       bx=0.75,
       by=lambda x: np.where(x < x0, 1.0, -1.0),
       p=lambda x: np.where(x < x0, 1.0, 0.1),
   )

   # 3. What to write to disk, and when
   ph.MHDDiagnostics(quantity="rho", write_timestamps=[final_time])

   Simulator(sim).run()

   # 4. Read the density back and plot it
   run = Run("mhd_run")
   run.GetMHDrho(final_time).plot(filename="rho.png")

Run it with ``python3 shock.py``. It takes a few seconds and writes
``rho.png``. You should see the density profile of the shock tube after
t = 40: flat plateaus separated by several waves that left the initial jump,
some going left and some going right. The box is periodic, so the box edges
are also a jump, and waves leave them too.

.. plot::
   :caption: The initial state of this example, drawn from the formulas (not
             from a PHARE run). The jump at the middle of the box is what
             produces the waves.

   import numpy as np
   import matplotlib.pyplot as plt

   plt.rcParams["font.size"] = 11
   L = 400.0
   x = np.linspace(0, L, 801)
   left = x < L / 2
   fig, axes = plt.subplots(3, 1, figsize=(3.4, 4.5), sharex=True)
   for ax, (name, l, r) in zip(
       axes, [("density", 1.0, 0.125), ("pressure", 1.0, 0.1), ("$B_y$", 1.0, -1.0)]
   ):
       ax.plot(x, np.where(left, l, r), color="k")
       ax.set_ylabel(name)
   axes[-1].set_xlabel("x")
   fig.tight_layout()


What it does
------------

An MHD run evolves one fluid (density, velocity, pressure) and the magnetic
field on a grid, with a :term:`finite volume` scheme. Each time step has three
choices you make in :py:class:`~pyphare.pharein.Simulation`:

.. plot::
   :include-source: False
   :caption: One step on four cells (made-up numbers): 1, each cell's limited slope gives the
             values on both sides of a face; 2, one flux per face; 3, each cell average changes
             by -(dt/dx) times its flux difference, here F3 - F2 for cell 2.

   from phare_draw import *
   u = [1.0, 0.9, 0.5, 0.4]            # cell averages (any conserved quantity)
   sl = [0.0, -0.1, -0.1, 0.0]         # limited slopes per cell (one cell = 1)
   Y = lambda v: 2.5 * (v - 0.3)         # value -> drawing height
   fig, (a, b, c) = strip(3, ["reconstruction: face values",
                              "Riemann solver: fluxes",
                              "time stepper: update"], height=1.9, aspect=False, phone=True)
   for ax in (a, b, c):
       for i in range(5):
           ax.plot([i, i], [0, 2.2], color=C["light"], lw=1)
       for i, v in enumerate(u):
           ax.plot([i, i + 1], [Y(v), Y(v)], color=C["grey"], lw=2)
           text(ax, i + 0.5, -0.3, f"{i}", color=C["grey"])
   for i, (v, s) in enumerate(zip(u, sl)):
       a.plot([i, i + 1], [Y(v - s / 2), Y(v + s / 2)], color=C["rho"], lw=2.2)
   mark(a, "primal", 2, Y(u[1] + sl[1] / 2), color=C["rho"], size=8)
   mark(a, "dual", 2, Y(u[2] - sl[2] / 2), color=C["rho"], size=8)
   text(a, 2.12, 2.1, "face values", ha="left", color=C["rho"])
   for f in (1, 2, 3):
       arrow(b, (f - 0.35, 2.05), (f + 0.35, 2.05), color=C["E"], lw=2)
       text(b, f, 2.35, f"F{f}", color=C["E"])
   for i in (1, 2):
       du = 0.04 if i == 1 else -0.04
       c.plot([i, i + 1], [Y(u[i] + du), Y(u[i] + du)], color=C["ok"], lw=2.4, ls="--")
   text(c, 2.5, 2.25, "F3 - F2", color=C["ok"])
   finish(fig, [a, b, c], max_h=1.7)

* :doc:`/theory/mhd`: the equations.
