=====================
Your first simulation
=====================

**Goal:** run a small hybrid simulation on your laptop (it takes well under a
minute), then plot one quantity and check that the physics is right.

The simulation is a single Alfvén wave travelling along a uniform magnetic
field, in one dimension. We know the answer in advance: the wave must move at
the Alfvén speed. That is what makes it a good first test.

.. contents:: On this page
   :local:
   :depth: 1


Before you start
----------------

You need PHARE built on your machine (see :doc:`install` and :doc:`build`), and
Python must be able to find two folders: PHARE's Python package ``pyphare``,
and your build folder (it holds the compiled part of PHARE).

.. code-block:: bash

   export PYTHONPATH=/path/to/PHARE/pyphare:/path/to/build:$PYTHONPATH

Replace ``/path/to/PHARE`` with the folder where you cloned PHARE, and
``/path/to/build`` with your build folder. You also need ``numpy``,
``matplotlib`` and ``h5py`` (installed with PHARE's Python requirements).

Work in an empty folder: the run writes its files where you start it.


Step 1 — write the input script
-------------------------------

A PHARE simulation is a Python script. Save this one as
``first_simulation.py``:

.. code-block:: python

   import numpy as np
   import pyphare.pharein as ph
   from pyphare.simulator.simulator import Simulator

   # 1. The simulation: grid, time steps, where to write the output
   sim = ph.Simulation(
       cells=200,                 # number of grid cells (one number = 1D)
       dl=0.5,                    # cell size, so the box is 200 * 0.5 = 100 long
       time_step=0.01,
       final_time=50.0,           # 5000 time steps
       diag_options={"format": "phareh5",
                     "options": {"dir": "first_run", "mode": "overwrite"}},
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
   ph.ElectromagDiagnostics(quantity="B", write_timestamps=[0.0, 25.0, 50.0])

   # 5. Run
   Simulator(sim).run()


Step 2 — run it
---------------

.. code-block:: bash

   python3 first_simulation.py

PHARE prints one line per time step, like ``t = 12.34000 - 0.00240sec - total
0:00:03.1``: the simulation time, the wall-clock time of that step, and the
total so far. The run ends with ``Finished at ...``. On a recent laptop this
takes 10 to 30 seconds.

Before the first time step you may also see a few messages, such as
``Warning: PHARE was not built with the configurator active``, or a note from
SAMRAI (the mesh library PHARE uses) starting with ``Allowing patches smaller
than the given smallest patch size``. They are harmless: they do not stop or
change the run.

When it is done, your folder contains:

.. code-block:: text

   first_simulation.py
   first_run/EM_B.h5      <- the magnetic field at t = 0, 25 and 50
   .log/                  <- messages from the compiled code, one file per process

``EM_B.h5`` is an `HDF5 <https://www.hdfgroup.org/solutions/hdf5/>`_ file.
You do not need to open it by hand: PHARE's analysis package reads it for you.


Step 3 — plot the wave
----------------------

Save this as ``plot_first_simulation.py`` in the same folder, and run it with
``python3 plot_first_simulation.py``:

.. code-block:: python

   import matplotlib.pyplot as plt
   from pyphare.pharesee.run import Run

   run = Run("first_run")                  # the folder written by the simulation

   fig, ax = plt.subplots()
   for t in run.times("B"):                # the 3 output times: 0, 25, 50
       B = run.GetB(t, merged=True)        # B at time t, on one grid
       by, (x,) = B["By"]                  # By as a function, and its x points
       ax.plot(x, by(x), label=f"t = {t:g}")

   ax.set_xlabel("x")
   ax.set_ylabel("By")
   ax.legend()
   fig.savefig("first_simulation.png")

``B["By"]`` gives two things: ``by``, a function that returns By at any
position, and a tuple with the grid coordinates, here just one array ``x``
since the run is 1D. ``by, (x,) = ...`` unpacks both in one line, and
``by(x)`` evaluates By on the grid.

Open ``first_simulation.png``. You should see three cosines, each shifted to
the right by about 25 compared to the previous one. Below is what the theory
predicts; your curves follow it, with small wiggles on top.

.. plot::
   :caption: Expected result: :math:`B_y(x) = 0.1\cos\big(k(x - t)\big)` at
             t = 0, 25 and 50. Drawn from the formula, not from a PHARE run.
             Your curves from PHARE have the same shape plus small wiggles.

   import numpy as np
   import matplotlib.pyplot as plt

   L, a = 100.0, 0.1
   k = 2 * np.pi / L
   x = np.linspace(0, L, 400)
   fig, ax = plt.subplots(figsize=(6, 3.2))
   for t in (0, 25, 50):
       ax.plot(x, a * np.cos(k * (x - t)), label=f"t = {t}")
       ax.axvline(t, color="0.7", lw=0.8, ls=":")
   ax.annotate("", xy=(50, 0.108), xytext=(0, 0.108),
               arrowprops=dict(arrowstyle="->", color="0.3"))
   ax.text(25, 0.111, "speed 1 (theory)", ha="center", va="bottom", fontsize=11)
   ax.set_ylim(-0.17, 0.135)
   ax.set_xlabel("x")
   ax.set_ylabel("By")
   ax.legend(loc="lower center", ncol=3, fontsize=11, frameon=False)
   fig.tight_layout()

The crest of the wave sits at x = 0, then about 25, then about 50: it moved
about 50 in a time of 50. Its speed is close to 1, the Alfvén speed in PHARE
units (next section). The particle noise makes the crest position uncertain
by about one cell either way.
**Your first simulation reproduces the right physics.**

The wiggles come from the particles: each cell holds 200 macro-particles (see
below), and a finite number of particles always gives some random noise. More
particles per cell means less noise and a longer run.


What just happened
------------------

The script declared four blocks, then handed them to the ``Simulator``, which
runs the compiled code and writes the output files. The analysis package,
``pharesee``, reads those files back. Click a box to read more.

.. mermaid::

   flowchart LR
     S["Simulation<br/>grid, time, output"] --> R
     M["MaxwellianFluidModel<br/>B field + protons"] --> R
     E["ElectronModel<br/>electron fluid"] --> R
     D["ElectromagDiagnostics<br/>what to write"] --> R
     R["Simulator(sim).run()"] --> F[("first_run/EM_B.h5")]
     F --> P["Run('first_run')<br/>pharesee"] --> G["plot"]
     click S href "../user_guide/simulation.html"
     click M href "../user_guide/hybrid_inputs.html"
     click E href "../user_guide/electrons.html"
     click D href "../user_guide/diagnostics.html"
     click R href "../user_guide/running.html"
     click P href "../data_analysis/reading_data.html"

**Hybrid model.** PHARE's default model is *hybrid*: ions are particles,
electrons are a fluid. That is why you gave a population of protons (with
``nbr_part_per_cell`` macro-particles in each cell, drawn from a Maxwellian
distribution of the given density, bulk velocity and thermal speed) and an
electron model. A macro-particle stands for many real ions.

**Units.** PHARE works in normalized units (details in
:doc:`/theory/normalization`). With a reference field :math:`B_0` and density
:math:`n_0`:

* lengths are in proton inertial lengths :math:`\delta_0` (also written
  :math:`d_i`), so the box is 100 :math:`d_i` long;
* times are in inverse proton gyrofrequencies :math:`1/\omega_0`;
* speeds are in Alfvén speeds :math:`v_0 = B_0/\sqrt{\mu_0 m_p n_0}`;
* temperatures (``Te`` here) are in :math:`m_p v_0^2`: the proton temperature
  is :math:`m_p v_{th}^2 = 0.01`, the same as ``Te``.

Here B = 1 and n = 1, so the Alfvén speed is exactly 1. The wavelength (100
:math:`d_i`) is much longer than :math:`d_i`, so the wave is very nearly a
plain Alfvén wave and moves at very nearly that speed (theory gives a few
percent less at this wavelength). For shorter wavelengths, a hybrid code shows the
dispersive behavior of the waves (see :doc:`/validation/hybrid_dispersion`).

**Hidden default.** PHARE adds a small hyper-resistivity (1e-4) by default. It
damps noise at the scale of one cell and leaves this long wave untouched; see
:doc:`/user_guide/electrons`.

**Why** ``vy = -by`` **and** ``vz = -bz``? An Alfvén wave going in the
:math:`+x` direction, along :math:`\mathbf{B}_0`, has its velocity perturbation
opposite to its magnetic perturbation (in these units). With the opposite sign,
the wave would go to the left.

**Periodic box.** With no boundary given, the box is periodic: what leaves on
the right comes back on the left. This is why the wavelength fits the box
exactly.

**merged=True.** PHARE stores its data in pieces of grid called *patches*
(see :doc:`/data_analysis/patch_hierarchy`). ``merged=True`` puts them back
together and gives you ``By`` as a function you can evaluate at the grid
points ``x``.


What each setting does
----------------------

.. list-table::
   :header-rows: 1
   :widths: 25 50 25

   * - Setting
     - Meaning
     - In this example
   * - ``cells``, ``dl``
     - number of cells and cell size; one number means 1D, a tuple of 2 or 3
       numbers means 2D or 3D. You can give ``domain_size`` instead of one of them.
     - 200 cells of 0.5
   * - ``time_step``, ``final_time``
     - time step and end time. You can give ``time_step_nbr`` instead of one of them.
     - 0.01 and 50
   * - ``diag_options``
     - output format and folder
     - ``first_run``
   * - ``bx``, ``by``, ``bz``
     - initial magnetic field, a number or a function of ``x``
     - uniform :math:`B_x` + wave
   * - ``protons={...}``
     - one ion population; the key (here ``protons``) is its name, you choose it
     - density 1, thermal speed 0.1
   * - ``nbr_part_per_cell``
     - macro-particles per cell for this population (default 100)
     - 200
   * - ``ElectronModel``
     - electron fluid. The *closure* is the rule that gives the electron
       pressure; ``isothermal`` (pressure = ``Te`` × density, ``Te`` fixed) is
       the only one. ``Te`` defaults to 0.1.
     - ``Te=0.01``
   * - ``ElectromagDiagnostics``
     - which field (``"E"`` or ``"B"``) to write, at which times
     - B at t = 0, 25, 50

All options are in the :doc:`/user_guide/simulation`,
:doc:`/user_guide/hybrid_inputs`, :doc:`/user_guide/electrons` and
:doc:`/user_guide/diagnostics` pages.


Try it yourself
---------------

Each change below takes one edit and one run.

* **Less noise:** set ``"nbr_part_per_cell": 1000``. The curves get smoother;
  the run takes about 5 times longer.
* **Wave going left:** change the signs of ``vy`` and ``vz``. The crest now
  moves to the left.
* **A shorter wave:** set ``L / 4`` as the wavelength (``k = 2 * np.pi / (L / 4)``).
  Four crests now fit in the box, and the wave is about 10% slower than
  speed 1: at t = 50 the crest that started at x = 0 is near x = 45, not 50.
  The wave speed now depends on the wavelength: this is dispersion.
* **Also write E:** add ``ph.ElectromagDiagnostics(quantity="E", write_timestamps=[0.0, 25.0, 50.0])``
  and read it with ``run.GetE(t, merged=True)``.


Common mistakes
---------------

* ``ModuleNotFoundError: No module named 'pyphare'`` — Python does not know
  where ``pyphare`` is. Add ``/path/to/PHARE/pyphare`` to ``PYTHONPATH``
  (see `Before you start`_).
* ``ModuleNotFoundError: No module named 'pybindlibs'`` (inside an error from
  ``Simulator.setup()``) — Python finds ``pyphare`` but not the compiled code.
  Add your build folder to ``PYTHONPATH``.
* ``RuntimeError: Error - no electrons registered to this Simulation`` — a
  hybrid simulation needs an ``ElectronModel`` block.
* ``Warning: Simulation is periodic but some functions are not`` — your initial
  profiles do not repeat exactly across the box (for example a wavelength that
  does not divide the box length). The run goes on, but the jump at the box
  edge can launch waves you did not ask for.
* Python ``from pyphare...`` works in one terminal but not in another — the
  ``export PYTHONPATH=...`` line only holds in the terminal where you typed it.
  Put it in your shell start-up file (for example ``~/.bashrc``) to keep it.


Next
----

:doc:`first_result` shows what is inside the output folder and how to plot
the other quantities PHARE can write.
