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
                     "options": {"dir": "first_run"}},   # no "mode": see the warning below
       restart_options={"dir": "restarts",
                        "mode": "conserve",
                        "restart_time": 25.0},           # read the restart at t = 25
   )

   # ... same model and electrons as run1.py ...

   ph.ElectromagDiagnostics(quantity="B", write_timestamps=[25.0, 50.0])

Run ``python3 run2.py``. It starts at t = 25, not at t = 0, and stops at
t = 50. Afterwards ``first_run/EM_B.h5`` holds the times 0, 25 and 50: t = 0
and 25 from the first run, and t = 50 from the second.

.. warning::

   Remove ``"mode": "overwrite"`` from ``diag_options`` in the script that
   restarts. If you keep it, the diagnostics files of the first run are emptied
   when the second run starts, without any message: ``EM_B.h5`` then holds only
   t = 25 and 50. Details in `Constraints`_.

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



How it works
------------

* **Writing.** At each time in ``timestamps`` PHARE writes a folder named after
  that time (``restarts/00025.00000``, five digits, a dot, five decimals).
  ``restart_options`` needs a ``mode``; without it PHARE stops.
* **Reading.** If ``restart_options`` contains ``restart_time``, PHARE loads
  the particles and fields saved at that time and goes on from there. Keep the model and electron
  blocks in the script: it is the same simulation, continued.
* **Check at start.** The new script must describe the same grid as the one
  that wrote the restart: ``cells``, ``dl`` and ``max_nbr_levels`` are compared,
  and PHARE stops if one differs. Nothing else is compared: keep every other
  setting identical too (model, populations, ``interp_order``, scheme options,
  ``time_step``); only ``final_time`` (or ``time_step_nbr``), the diagnostics
  and the restart options should change.
* Give ``timestamps`` as multiples of ``time_step`` (the simulation only
  passes through multiples of ``time_step``), in ascending order, between the
  start time and ``final_time``.


Options
-------

All keys go in the ``restart_options`` dictionary of
:class:`~pyphare.pharein.Simulation`. Any other key is rejected.

.. list-table::
   :header-rows: 1
   :widths: 20 60 20

   * - Option
     - Meaning
     - Default
   * - ``mode``
     - ``"conserve"`` or ``"overwrite"``. With ``"conserve"``, a restart time
       for which a folder already exists is skipped (not written again). With
       ``"overwrite"``, existing restart files are replaced.
     - none: **required**
   * - ``dir``
     - Folder holding the restart folders, for writing and for reading.
       **Always set it** when you use ``"auto"`` or ``keep_last``: without
       ``dir``, those two look in the current folder, not in ``phare_outputs``,
       so ``"auto"`` never finds a restart and ``keep_last`` deletes nothing.
     - ``"phare_outputs"``
   * - ``timestamps``
     - Simulation times at which to write a restart. A list or a numpy array.
     - none: no restart is written
   * - ``elapsed_timestamps``
     - List of *wall-clock* times since the run started (seconds or
       ``datetime.timedelta``), ascending; one restart at each.
       Can be combined with ``timestamps``.
     - none
   * - ``restart_time``
     - Simulation time of the restart to read, or ``"auto"`` for the latest
       one in ``dir``. If ``"auto"`` finds nothing, the run starts from
       scratch.
     - none: normal start
   * - ``keep_last``
     - Keep only the ``N`` most recent restart folders in ``dir``; older ones
       are deleted as new ones are written.
     - none: keep all

.. tip::

   On a cluster you rarely know which simulation time the job will have
   reached when it is killed. Prefer ``elapsed_timestamps`` (for example 10 h
   and 23 h, as ``datetime.timedelta(hours=10)``) so that a restart is written
   in time, whatever the speed of the run. Folders are still named after the
   *simulation* time at which the restart was written (for instance
   ``00011.46000``).

.. tip::

   One script for every job: give ``restart_time="auto"``, a ``dir``, and a
   ``timestamps`` or ``elapsed_timestamps`` list. The first job finds no restart and starts
   from scratch; each later job continues from the latest one. Add
   ``keep_last=2`` to limit the disk used.


Constraints
-----------

**Same number of MPI ranks.** The restart is written as one file per rank,
inside a folder named after the number of ranks (``nodes.0000001`` for one
rank, ``nodes.0000002`` for two, see ``restarts/00025.00000/restore.000000/``).
A restart written with N ranks can only be read with N ranks.

**Restart folder and diagnostics folder have different default names.** If you
give no ``dir``, restarts go to ``phare_outputs`` and diagnostics go to
``phare_output`` (no ``s``) when no ``options`` dictionary is given in
``diag_options``. Setting both ``dir`` entries explicitly avoids any confusion.

**Diagnostics files when restarting.** The ``mode`` of ``diag_options["options"]``
only accepts ``"overwrite"``. If you leave it out, existing diagnostics files
are opened and the new times are added at the end, which is what you want in a
restart. If you keep ``"overwrite"`` in the second script, the files written by
the first run are wiped at the start: for example, ``EM_B.h5`` then held only
t = 25 and 50, the t = 0 snapshot was lost. Remove ``"mode": "overwrite"`` (or
use another ``dir``) in the script that restarts.


Common mistakes
---------------

* ``ValueError: Restart mode not set, valid modes are ['conserve', 'overwrite']``
  — ``restart_options`` has no ``mode``. Add ``"mode": "conserve"``.
* ``ValueError: invalid option (mod), valid options are ['dir', 'elapsed_timestamps', 'timestamps', 'mode', 'restart_time', 'keep_last']``
  — a typo in a key name (here ``mod`` for ``mode``).
* ``ValueError: PHARE restart file not found for time 30.0`` — there is no
  folder ``restarts/00030.00000``. Look at what is in ``dir``, or use
  ``"restart_time": "auto"``.
* ``ERROR: Simulation not restartable - variable mismatch`` then
  ``cells current([100]) ne previous([200])`` and
  ``ValueError: deserialized Restart simulation is incompatible with configured simulation parameters``
  — the script changed ``cells`` (or ``dl`` or ``max_nbr_levels``). Use the
  values of the run that wrote the restart.
* ``RuntimeError: Unable to open file restarts/00025.00000/restore.000000/nodes.0000002/proc.0000000 (File accessibility) Unable to open file``
  (together with an HDF5 ``No such file or directory`` message) — you started
  with ``mpirun -n 2`` a restart written with one rank. Use the same number of
  ranks as in the first run.


Next
----

:doc:`running` explains how to launch a simulation, including with MPI. How
this works inside: :doc:`/dev/diagnostics_restarts`.
