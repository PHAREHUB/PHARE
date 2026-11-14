========
Restarts
========

**Goal:** stop a simulation, then continue it later from where it stopped.

On a cluster, jobs are killed after a fixed wall-clock time (for example 24
hours), and a long simulation rarely fits in one job. The answer is to write
**restart files** during the run (a full snapshot of the particles and fields),
and to start the next job from one of them. A restart is not the same as the
diagnostics: diagnostics are what you analyse (see :doc:`diagnostics`), restart
files are only there so PHARE can resume.

.. contents:: On this page
   :local:
   :depth: 1


Minimal example: stop at t = 25, continue to t = 50
---------------------------------------------------

Two scripts. The first one is the :doc:`/getting_started/first_simulation`
script, shortened to ``final_time=25`` and with a ``restart_options`` that
writes one restart at t = 25. Save it as ``run1.py``:

.. code-block:: python

   import numpy as np
   import pyphare.pharein as ph
   from pyphare.simulator.simulator import Simulator

   # 1. The simulation: grid, time steps, where to write the output
   sim = ph.Simulation(
       cells=200,                 # number of grid cells (one number = 1D)
       dl=0.5,                    # cell size, so the box is 200 * 0.5 = 100 long
       time_step=0.01,
       final_time=25.0,           # 2500 time steps
       diag_options={"format": "phareh5",
                     "options": {"dir": "first_run", "mode": "overwrite"}},
       restart_options={"dir": "restarts",
                        "mode": "conserve",
                        "timestamps": [25.0]},     # write a restart at t = 25
   )

   # 2. The initial plasma: uniform protons and a field B0 = 1 along x,
   #    plus a small circularly polarized Alfven wave
   L = 100.0                      # box length
   k = 2 * np.pi / L              # one wavelength in the box
   a = 0.1                        # wave amplitude

   def density(x): return 1.0
   def bx(x):      return 1.0
   def by(x):      return a * np.cos(k * x)
   def bz(x):      return a * np.sin(k * x)
   def vx(x):      return 0.0
   def vy(x):      return -a * np.cos(k * x)
   def vz(x):      return -a * np.sin(k * x)
   def vth(x):     return 0.1     # proton thermal speed

   ph.MaxwellianFluidModel(
       bx=bx, by=by, bz=bz,
       protons={"charge": 1, "density": density,
                "vbulkx": vx, "vbulky": vy, "vbulkz": vz,
                "vthx": vth, "vthy": vth, "vthz": vth,
                "nbr_part_per_cell": 200},
   )

   # 3. The electrons: a fluid with a fixed temperature
   ph.ElectronModel(closure="isothermal", Te=0.01)

   # 4. What to write to disk, and when
   ph.ElectromagDiagnostics(quantity="B", write_timestamps=[0.0, 25.0])

   # 5. Run
   Simulator(sim).run()

Run it with ``python3 run1.py`` (about 7 seconds). It creates two folders:
``first_run/`` (diagnostics) and ``restarts/00025.00000/`` (the restart).

The second script, ``run2.py``, is the same file with three changes. The grid,
the plasma and the electrons must stay exactly as they were: it is the same
simulation, continued.

.. code-block:: python

   sim = ph.Simulation(
       cells=200,
       dl=0.5,
       time_step=0.01,
       final_time=50.0,           # was 25.0: now run up to t = 50
       diag_options={"format": "phareh5",
                     "options": {"dir": "first_run"}},   # no "mode": see below
       restart_options={"dir": "restarts",
                        "mode": "conserve",
                        "restart_time": 25.0},           # read the restart at t = 25
   )

   # ... same model and electrons as run1.py ...

   ph.ElectromagDiagnostics(quantity="B", write_timestamps=[25.0, 50.0])

Run ``python3 run2.py``. It starts at t = 25, not at t = 0, and stops at
t = 50. Afterwards ``first_run/EM_B.h5`` holds the times 0, 25 and 50: t = 0
and 25 from the first run, and t = 50 from the second.

.. plot::
   :include-source: False
   :caption: Run 1 (blue) writes B at t = 0 and 25 (squares) and the restart files at 25;
             run 2 (orange) starts from them and adds B at t = 50 to the same file.

   from phare_draw import *
   from matplotlib.patches import Rectangle
   fig, (a,) = strip(1, ["Two jobs, one simulation"], height=2.6, phone=True)
   s = 0.2  # 1 drawing unit = 5 time units
   def bar(x0, x1, y, color, label):
       a.add_patch(Rectangle((x0*s, y), (x1-x0)*s, 0.5, fc=color, ec=color, alpha=0.25, lw=1.5))
       text(a, (x0+x1)/2*s, y+0.25, label, color=C["ink"])
   arrow(a, (0, 0), (10.8, 0), color=C["grey"], lw=1.2)
   text(a, 11.2, 0, "t", color=C["grey"])
   for t in (0, 25, 50):
       a.plot([t*s, t*s], [-0.12, 0.12], color=C["grey"], lw=1.2)
       text(a, t*s, -0.45, f"{t}", color=C["grey"])
   bar(0, 25, 1.2, C["L1"], "run1.py")
   bar(25, 50, 2.6, C["L2"], "run2.py")
   mark(a, "primal", 25*s, 2.05, color=C["ink"], size=10)
   text(a, 25*s + 0.35, 2.05, "restart", ha="left")
   arrow(a, (25*s, 2.2), (25*s, 2.55), color=C["ink"], lw=1.5)
   for t, c in ((0, C["L1"]), (25, C["L1"]), (50, C["L2"])):
       mark(a, "dual", t*s, 0.55, color=c, size=9)
   finish(fig, [a])

this works inside: :doc:`/dev/diagnostics_restarts`.
