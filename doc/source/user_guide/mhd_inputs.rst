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
       domain_boundaries={"xlower": {"type": "open"},   # outflow on both sides
                          "xupper": {"type": "open"}},  # (the default is periodic)
       diag_options={"format": "phareh5",
                     "options": {"dir": "mhd_run", "mode": "overwrite"}},
   )

   # 2. The initial fluid: left and right states
   ph.MHDModel(
       density=lambda x: np.where(x < x0, 1.0, 0.125),
       vx=0.0,                     # the default is 1, not 0: see below
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
some going left and some going right. The two ``open`` sides let them leave
the box. Without ``domain_boundaries`` the box would be periodic: its edges would
then be a second jump, with waves leaving it too.

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


Scheme options
--------------

The numerical scheme is chosen by these arguments of ``Simulation`` (details and
consequences: :doc:`/numerics/mhd_finite_volume`):

* ``reconstruction``: ``"Constant"``, ``"Linear"``, ``"WENO3"``, ``"WENOZ"`` or ``"MP5"``
  (:doc:`/numerics/mhd_reconstruction`);
* ``limiter``: ``"None"``, ``"VanLeer"`` or ``"MinMod"``; only ``"Linear"`` uses one;
* ``riemann``: ``"Rusanov"``, ``"HLL"`` or ``"HLLD"`` (:doc:`/numerics/mhd_riemann_solvers`);
* ``mhd_timestepper``: ``"Euler"``, ``"TVDRK2"``, ``"TVDRK3"`` or ``"SSPRK4_5"``, mandatory
  (:doc:`/numerics/mhd_time_integration`);
* ``hall``: ``True`` or ``False`` (default ``False``);
* ``gamma``: adiabatic index, greater than 1 (default 5/3);
* ``eta``, ``nu``: resistivity and hyper-resistivity (default 0; ``nu`` needs ``hall=True``),
  with ``hyper_mode``: ``"constant"`` (default) or ``"spatial"``.

``reconstruction``, ``limiter``, ``riemann``, ``hall`` and the dimension are fixed when PHARE
is compiled: not every combination is built by default (:doc:`/numerics/mhd_finite_volume`).

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

* the **reconstruction** (with its **limiter**) estimates the fluid on the two
  sides of every cell face;
* the **Riemann solver** turns those two states into a flux through the face;
* the **time stepper** combines the fluxes into a new state.

How these work is described in :doc:`/numerics/mhd_finite_volume`.

Units
^^^^^

MHD uses the same :doc:`/theory/normalization` as the hybrid model. Lengths are
in the proton inertial length, times in the inverse proton gyro-frequency,
the magnetic field in a reference field :math:`B_0`, the mass density in
:math:`m_p n_0`, and the pressure in :math:`B_0^2/\mu_0`. ``p`` is the thermal
pressure of the fluid only (the magnetic pressure is not included in it), and
``gamma`` links it to the internal energy: :math:`e = p/(\gamma - 1)`.


Describing the initial fluid
----------------------------

:py:class:`~pyphare.pharein.MHDModel` takes one profile per quantity. Each is
either a number (uniform) or a function of the coordinates: ``f(x)`` in 1D,
``f(x, y)`` in 2D, ``f(x, y, z)`` in 3D. The function must return an array of the
same size as its input (use ``np.where``, or ``0 * x + value``).

.. list-table::
   :header-rows: 1

   * - Option
     - Meaning
     - Default
   * - ``density``
     - mass density
     - 1
   * - ``vx``, ``vy``, ``vz``
     - velocity components
     - 1, 0, 0
   * - ``bx``, ``by``, ``bz``
     - magnetic field components
     - 1, 0, 0
   * - ``p``
     - thermal pressure
     - 1

.. warning::

   The defaults are not at rest: ``vx`` is 1 and ``bx`` is 1 unless you say
   otherwise. This is why the example above gives ``vx=0.0``: without that
   line the whole shock tube is carried along x at speed 1, with no message.

Only one model can be created per simulation, and it must come after
:py:class:`~pyphare.pharein.Simulation`.

The magnetic field you give should have no divergence (:math:`\nabla\cdot B = 0`);
a uniform component, as above, is always safe.


More on the scheme options
--------------------------

* ``reconstruction``, ``limiter`` and ``riemann`` have no default (an empty
  name) and are **case-sensitive**: write them exactly as in the list above.
  ``mhd_timestepper`` accepts any letter case.
* In an MHD run the hyper-resistivity is ``nu`` and the resistivity ``eta``. The
  hybrid options ``resistivity`` and ``hyper_resistivity`` are accepted in an
  MHD run and have no effect.
* ``max_mhd_level`` is meant for a future run mixing MHD and hybrid
  :term:`levels <level>`, which does not work yet: leave it at its default.
* Which combinations of ``reconstruction``, ``limiter``, ``riemann``, ``hall``
  and dimension are built by default: :ref:`which-combinations`. To see what
  your own build contains, run
  ``python3 -c "from pyphare import cpp; print(cpp.built_simulator_ids())"``
  (MHD entries read ``<dim>_<reconstruction>_<limiter>_<riemann>_<hall>``).


Domain boundaries
-----------------

By default the box is periodic. Walls, outflow and inflow are set with
``domain_boundaries``, which only MHD supports: see
:doc:`/user_guide/boundary_conditions`. The example above uses two ``open``
(outflow) sides.


Diagnostics
-----------

``MHDDiagnostics`` writes ``rho``, ``V``, ``P``, ``rhoV`` and ``Etot``; the
magnetic field is written with ``ElectromagDiagnostics(quantity="B", ...)``
(see :doc:`diagnostics`). The electric field is **not** written in an MHD run:
``ElectromagDiagnostics(quantity="E", ...)`` stops the run at the start (see
below).

They are read with ``run.GetMHDrho``, ``GetMHDV``, ``GetMHDP``, ``GetMHDrhoV``,
``GetMHDEtot`` and ``run.GetB`` (see :doc:`/data_analysis/reading_data`).


Common mistakes
---------------

* **Combination not built** (a typo, wrong letter case, or a valid scheme that
  was not compiled). Wrong case, for instance ``reconstruction="linear"``,
  gives:

  .. code-block:: text

     ModuleNotFoundError: No module named 'pybindlibs.cpp_1_linear_VanLeer_Rusanov_false'.
     This run needs the build permutation
         1,linear,VanLeer,Rusanov,false
     Add it to a res/sim permutation file (res/sim/all.txt) and rebuild.
     Permutations currently built:
         ...

  Fix: use the exact names above, and a combination that is built. The Python
  side does not check the names, so the same message appears for
  ``riemann="rusanov"``, or for a limiter with a non-Linear reconstruction such
  as ``reconstruction="WENOZ", limiter="VanLeer"``
  (``1,WENOZ,VanLeer,Rusanov,false``). Leaving ``limiter`` out gives the same
  message with an empty name: write ``limiter="None"`` for the schemes that
  take no limiter.

* **No** ``ph.MHDModel(...)`` **block**: ``AttributeError: 'NoneType' object has
  no attribute 'model_dict'`` from ``Simulator.setup()``. Add the block after
  ``ph.Simulation(...)``.

* **Forgetting** ``vx`` (or ``bx``): the defaults are 1, not 0 (see the warning
  above). The run works, with a flow you did not ask for.

* **Asking for E**:

  .. code-block:: text

     RuntimeError: E not currently available in MHD diagnostics

  followed by ``Error in Simulator.initialize()``. Fix: in an MHD run write
  only ``ElectromagDiagnostics(quantity="B", ...)``.

* ``nu`` **without the Hall term**:

  .. code-block:: text

     RuntimeError: Error - hyper-resistivity (nu > 0) requires the Hall term to be enabled

  Fix: add ``hall=True`` (and a combination built with ``hall`` true, such as
  1D ``Linear``/``VanLeer``/``Rusanov``), or set ``nu=0``.

* ``resistivity`` **or** ``hyper_resistivity`` **in an MHD run.** Both are hybrid
  options: an MHD run accepts them and ignores them. Use ``eta`` and ``nu``.

* **No time stepper**:

  .. code-block:: text

     ValueError: Error: MHDModel requires a non-empty mhd_timestepper

* **Unknown time stepper**, for example ``"rk4"``:

  .. code-block:: text

     ValueError: TimeIntegratorType: unknown value 'rk4', expected one of ['euler', 'tvdrk2', 'tvdrk3', 'ssprk4_5']

* ``gamma`` **of 1 or less**:

  .. code-block:: text

     ValueError: 'gamma' must be a scalar greater than 1, got 1.0

* ``max_mhd_level`` **too large**:

  .. code-block:: text

     ValueError: Error: max_mhd_level(3) should be less or equal to max_nbr_levels(1)

* **Using the hybrid model in an MHD run.** ``MaxwellianFluidModel`` belongs to the
  hybrid model. In a run declared with ``model_options=["MHDModel"]`` it fails
  with an unhelpful message as soon as it is declared:

  .. code-block:: text

     ValueError: HybridGridLayoutFor: interp_order is 0 (no hybrid model) and no explicit ghosts_nbr was given

  Fix: use :py:class:`~pyphare.pharein.MHDModel`.

* **Hybrid options in an MHD run.** ``interp_order`` is accepted by Python, but
  it asks for a hybrid+MHD module that is not built: the run stops with
  ``No module named 'pybindlibs.cpp_1_1_2_Linear_VanLeer_Rusanov_false'``.
  ``refined_particle_nbr`` alone is rejected at once: ``Error:
  refined_particle_nbr is meaningless without a hybrid model (no interp_order
  given), it should be omitted``. Leave both out.


Going further
-------------

* :doc:`/user_guide/boundary_conditions`: open, reflective and inflow boundaries.
* :doc:`/user_guide/diagnostics`: choosing what to write and when.
* :doc:`/numerics/mhd_finite_volume`: how the scheme works.
* :doc:`/theory/mhd`: the equations.


Reference
---------

.. autoclass:: pyphare.pharein.MHDModel
