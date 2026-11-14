=============================
Hybrid inputs: ions and field
=============================

**Goal:** describe the initial state of a hybrid run: the magnetic field and
one or more ion populations, each loaded as macro-particles.

In a hybrid run, ions are particles and electrons are a fluid. The ions and the
initial magnetic field are set with one block, ``MaxwellianFluidModel``; the
electrons with ``ElectronModel`` (see :doc:`electrons`). Declare
:doc:`Simulation <simulation>` first: the model needs to know the dimension.
Hybrid and MHD cannot be combined in one run yet (see :doc:`choosing_a_model`).

.. contents:: On this page
   :local:
   :depth: 1


Minimal example: protons and alpha particles
--------------------------------------------

A 1D box with protons that have a density bump in the middle, and 5 % of alpha
particles (charge 2, mass 4) drifting along the field. It runs in a few
seconds.

.. code-block:: python

   import numpy as np
   import pyphare.pharein as ph
   from pyphare.simulator.simulator import Simulator

   sim = ph.Simulation(
       cells=200, dl=0.5,
       time_step=0.01, final_time=1.0,
       diag_options={"format": "phareh5",
                     "options": {"dir": "two_species", "mode": "overwrite"}},
   )

   L = 100.0

   def n_protons(x):                    # a density bump in the middle of the box
       return 1.0 + 0.5 * np.exp(-((x - L / 2) / 5.0) ** 2)

   def n_alphas(x):
       return 0.05                      # 5 % of alpha particles, uniform

   ph.MaxwellianFluidModel(
       bx=1.0,                          # a number: uniform field
       protons={"charge": 1, "mass": 1, "density": n_protons,
                "vthx": 0.1, "vthy": 0.1, "vthz": 0.1,
                "nbr_part_per_cell": 100},
       alphas={"charge": 2, "mass": 4, "density": n_alphas,
               "vbulkx": 0.2,           # alphas drift along B
               "vthx": 0.05, "vthy": 0.05, "vthz": 0.05,
               "nbr_part_per_cell": 50},
   )
   ph.ElectronModel(closure="isothermal", Te=0.05)

   times = [0.0, 1.0]
   ph.FluidDiagnostics(quantity="density", population_name="protons", write_timestamps=times)
   ph.FluidDiagnostics(quantity="density", population_name="alphas", write_timestamps=times)
   ph.FluidDiagnostics(quantity="charge_density", write_timestamps=times)

   Simulator(sim).run()

Every keyword of ``MaxwellianFluidModel`` other than ``bx``, ``by``, ``bz`` is a
population: its name is yours to choose (``protons``, ``alphas``, ``beam``...)
and its value is a dictionary of parameters. Anything you leave out takes its
default value (table below).

You can check that PHARE adds the populations as expected: the ion charge
density is the sum of the population densities weighted by their charge,
:math:`\rho_c = n_p + 2 n_\alpha`. First copy the ``profile`` function (and
its imports) from :doc:`/getting_started/first_result`, then:

.. code-block:: python

   run = Run("two_species")
   x, n_p = profile(run.GetN(0.0, "protons"), "value")
   x, n_a = profile(run.GetN(0.0, "alphas"), "value")
   x, rho_c = profile(run.GetNi(0.0), "value")
   print(abs(rho_c - (n_p + 2 * n_a)).max())    # about 1e-7: single-precision output


How particles are loaded
------------------------

PHARE evaluates your profile functions at the center of every cell. In each
cell, it then creates ``nbr_part_per_cell`` macro-particles for each
population:

* their **positions** are drawn at random, uniformly inside the cell;
* their **velocities** are drawn from a Maxwellian: in each direction, a
  Gaussian centered on the bulk velocity (``vbulkx``...) with standard
  deviation the thermal velocity (``vthx``...);
* their **weight** is the density divided by the number of particles per cell.
  A cell where the density is twice as high has the same number of particles,
  each twice as heavy.

.. plot::
   :caption: Sketch of the loading (random numbers, not a PHARE run). Top: a
             density profile, evaluated at cell centers (dots). Bottom: each
             cell holds the same number of particles (here 4), placed at random
             in the cell; the marker area is the weight, density / 4.

   import numpy as np
   import matplotlib.pyplot as plt

   plt.rcParams["font.size"] = 11
   rng = np.random.default_rng(1)
   ncells, ppc, dx = 12, 4, 1.0
   xc = (np.arange(ncells) + 0.5) * dx
   n = 1.0 + 0.8 * np.exp(-((xc - 6.0) / 2.0) ** 2)
   fig, (a1, a2) = plt.subplots(nrows=2, sharex=True, figsize=(6, 4))
   xs = np.linspace(0, ncells * dx, 300)
   a1.plot(xs, 1.0 + 0.8 * np.exp(-((xs - 6.0) / 2.0) ** 2), color="0.6")
   a1.plot(xc, n, "o", color="C0")
   a1.set_ylabel("density")
   for i in range(ncells):
       a2.axvline(i * dx, color="0.85", lw=0.8)
       xp = i * dx + rng.uniform(0, dx, ppc)
       yp = rng.uniform(0.1, 0.9, ppc)
       a2.scatter(xp, yp, s=60 * n[i] / ppc * 4 / 2, color="C1")
   a2.set_yticks([])
   a2.set_xlabel("x (cells)")
   a2.set_ylabel("particles")
   fig.tight_layout()

Consequences worth knowing:

* A population with a low density (here the alphas) still has
  ``nbr_part_per_cell`` particles per cell, each lighter. You do not lose
  resolution on minor species.
* The noise of a moment (density, velocity...) goes down as one over the
  square root of the number of particles per cell: four times more particles,
  half the noise, about four times the cost.
* Cells where the density is below ``density_cut_off`` get no particles.
  Keep it far below the lowest density of your run: a hybrid run cannot have
  empty cells, because Ohm's law divides by the density (a run with an empty
  region crashes at the start with a segmentation fault). For a low-density region, use
  a small but non-zero density.
* Without a seed, each run draws different particles. With
  ``"init": {"seed": 42}``, the same input gives the same particles whatever
  the number of :term:`MPI` ranks or the patch layout (for a given compiler
  and standard library).


Profiles: numbers or functions
------------------------------

``bx``, ``by``, ``bz`` and every population parameter except ``charge``,
``mass``, ``nbr_part_per_cell``, ``density_cut_off`` and ``init`` can be:

* a **number**, for a uniform value (``bx=1.0``);
* a **function** of the coordinates, with one argument per dimension:
  ``f(x)`` in 1D, ``f(x, y)`` in 2D, ``f(x, y, z)`` in 3D. PHARE calls it with
  numpy arrays, so use ``np.cos``, ``np.exp``... rather than ``math.cos``.

Functions can use any Python you like, including helpers and values from the
``Simulation`` object. For example, two current sheets at :math:`y = 0.3 L_y`
and :math:`y = 0.7 L_y`, each of half-width 0.5, in 2D (a fragment: it goes in
a 2D script, after ``sim = ph.Simulation(...)``):

.. code-block:: python

   import numpy as np

   def density(x, y):
       Ly = sim.simulation_domain()[1]
       return (
           0.2
           + 1.0 / np.cosh((y - Ly * 0.3) / 0.5) ** 2
           + 1.0 / np.cosh((y - Ly * 0.7) / 0.5) ** 2
       )

   def S(y, y0, l):
       return 0.5 * (1.0 + np.tanh((y - y0) / l))

   def bx(x, y):
       Ly = sim.simulation_domain()[1]
       v1, v2 = -1.0, 1.0
       return v1 + (v2 - v1) * (S(y, Ly * 0.3, 0.5) - S(y, Ly * 0.7, 0.5))

A complete input built on these profiles (with a different background
density and a small perturbation) is ``tests/functional/harris/harris_2d.py``
(see :doc:`/examples/harris`).

**Periodic box.** Unless you set boundaries (only possible for MHD, see
:doc:`boundary_conditions`), the box is periodic, and every profile must take
the same value on both sides of the box. PHARE checks it when the model is
declared: if a profile does not repeat, it prints
``Warning: Simulation is periodic but some functions are not`` and goes on
(with ``Simulation(strict=True)``, it stops instead). The jump at the edge would
otherwise launch waves you did not ask for.

**Magnetic field and div B.** PHARE does not check that your initial field has
zero divergence: choose ``bx``, ``by``, ``bz`` so that
:math:`\nabla\cdot\mathbf{B} = 0`.


Options
-------

Field (keywords of ``MaxwellianFluidModel``):

.. list-table::
   :header-rows: 1
   :widths: 25 55 20

   * - Option
     - Meaning
     - Default
   * - ``bx``, ``by``, ``bz``
     - initial magnetic field, number or function
     - ``bx=1``, ``by=bz=0``

Population (keys of each population dictionary):

.. list-table::
   :header-rows: 1
   :widths: 25 55 20

   * - Key
     - Meaning
     - Default
   * - ``charge``
     - charge of one ion, in units of the proton charge
     - 1
   * - ``mass``
     - mass of one ion, in proton masses
     - 1
   * - ``density``
     - number density, number or function
     - 1
   * - ``vbulkx``, ``vbulky``, ``vbulkz``
     - bulk velocity, number or function
     - 0
   * - ``vthx``, ``vthy``, ``vthz``
     - thermal velocity: standard deviation of the Maxwellian in each
       direction, number or function
     - **1**
   * - ``nbr_part_per_cell``
     - macro-particles per cell
     - 100
   * - ``density_cut_off``
     - cells with a density below this get no particles
     - 1e-16
   * - ``init``
     - ``{"seed": <int>}`` to make the particle draw reproducible
     - random

Units are those of :doc:`/theory/normalization`: densities in :math:`n_0`,
velocities in Alfvén speed :math:`v_0`, magnetic field in :math:`B_0`. The
thermal velocity relates to the temperature by :math:`v_{th}^2 = T/m`, in these
units.


Common mistakes
---------------

* **Forgetting** ``vthx``, ``vthy``, ``vthz``. The default thermal velocity is
  **1**, the Alfvén speed: a hot plasma (:math:`\beta_i = 2` for density 1 and
  B = 1). Always give them.
* ``RuntimeError: A simulation must be declared before a model`` — declare
  ``ph.Simulation(...)`` before ``MaxwellianFluidModel``.
* ``RuntimeError: A model is already created`` — there can be only one
  ``MaxwellianFluidModel`` in a script. Put all populations in the same call.
* ``Warning: Simulation is periodic but some functions are not`` — a profile
  does not repeat across the box (see above).
* A segmentation fault (``Caught signal 11``) at the start of the run — part
  of the box has no particles (density zero, or below ``density_cut_off``).
  Give every cell a non-zero density.
* ``TypeError: MaxwellianFluidModel.add_population() got an unexpected keyword
  argument 'densty'. Did you mean 'density'?`` — a typo in a population key.
* A function with the wrong number of arguments (``f(x)`` in a 2D run) stops
  with an ``AssertionError`` when the model is declared.
* A run that stops with ``RuntimeError: Error - no electrons registered to this
  Simulation`` — add an :doc:`ElectronModel <electrons>`.


Reference
---------

.. autoclass:: pyphare.pharein.MaxwellianFluidModel
