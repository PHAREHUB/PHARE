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
   fig, axes = plt.subplots(3, 1, figsize=(6, 5.5), sharex=True)
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

.. mermaid::

   flowchart LR
     A["reconstruction (+ limiter)<br/>values on each cell face"] --> B["riemann<br/>flux through each face"]
     B --> C["mhd_timestepper<br/>advance in time"]

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

Note that the defaults are not at rest: ``vx`` is 1 and ``bx`` is 1 unless you
say otherwise. Only one model can be created per simulation, and it must come
after :py:class:`~pyphare.pharein.Simulation`.

The magnetic field you give should have no divergence (:math:`\nabla\cdot B = 0`);
a uniform component, as above, is always safe.


Simulation options for MHD
--------------------------

.. list-table::
   :header-rows: 1

   * - Option
     - Meaning
     - Default
   * - ``model_options``
     - ``["MHDModel"]`` for an MHD run
     - ``"HybridModel"``
   * - ``mhd_timestepper``
     - time integration scheme: ``"Euler"``, ``"TVDRK2"``, ``"TVDRK3"`` or ``"SSPRK4_5"`` (any letter case)
     - none: required
   * - ``reconstruction``
     - ``"Constant"``, ``"Linear"``, ``"WENO3"``, ``"WENOZ"`` or ``"MP5"`` (case-sensitive)
     - ``""``: you must set it
   * - ``limiter``
     - ``"VanLeer"`` or ``"MinMod"`` (only with ``"Linear"``), otherwise ``"None"`` (case-sensitive)
     - ``""``: you must set it
   * - ``riemann``
     - ``"Rusanov"``, ``"HLL"`` or ``"HLLD"`` (case-sensitive)
     - ``""``: you must set it
   * - ``hall``
     - include the Hall term in Ohm's law
     - ``False``
   * - ``gamma``
     - adiabatic index, must be larger than 1
     - 5/3
   * - ``eta``
     - resistivity
     - 0
   * - ``nu``
     - hyper-resistivity (needs ``hall=True``)
     - 0
   * - ``hyper_mode``
     - how ``nu`` is applied: ``"constant"`` or ``"spatial"``
     - ``"constant"``
   * - ``max_mhd_level``
     - meant for a future run mixing MHD and hybrid :term:`levels <level>`,
       which does not work yet: leave it at its default (at most ``max_nbr_levels``)
     - ``max_nbr_levels`` in a pure MHD run

In an MHD run the hyper-resistivity is ``nu``. The option ``hyper_resistivity``
belongs to the hybrid model: if you pass it in an MHD run it is accepted and
has no effect.

The scheme names are explained in :doc:`/numerics/mhd_finite_volume`.

Which combinations are available
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``reconstruction``, ``limiter``, ``riemann``, the dimension and ``hall`` are
compiled into PHARE: each combination is a separate module, and only some are
built by default. The ones below are built (``res/sim/all.txt``); anything else
needs a rebuild with an extra line in that file. ``Linear`` is the only
reconstruction that takes a limiter. In every row the time stepper is free to
choose, since it is not compiled in.

.. list-table::
   :header-rows: 1

   * - Dimension
     - reconstruction
     - limiter
     - riemann
     - hall
   * - 1D
     - Linear
     - VanLeer
     - Rusanov
     - False and True
   * - 1D
     - Linear
     - MinMod
     - Rusanov
     - False
   * - 1D
     - Constant
     - None
     - Rusanov
     - False
   * - 1D
     - Constant
     - None
     - HLL, HLLD
     - False
   * - 1D
     - WENO3, WENOZ, MP5
     - None
     - Rusanov
     - False
   * - 2D
     - Linear
     - VanLeer
     - Rusanov
     - False
   * - 2D
     - WENOZ
     - None
     - Rusanov
     - False and True
   * - 3D
     - Linear
     - VanLeer
     - Rusanov
     - False
   * - 3D
     - WENO3, MP5
     - None
     - Rusanov
     - False
   * - 3D
     - WENOZ
     - None
     - Rusanov
     - False and True

``HLL`` and ``HLLD`` are only built for 1D ``Constant``. To see what your
build contains, run
``python3 -c "from pyphare import cpp; print(cpp.built_simulator_ids())"``
(MHD entries read ``<dim>_<reconstruction>_<limiter>_<riemann>_<hall>``), or
simply try your combination: the error message lists what is built (see below).


Domain boundaries
-----------------

By default the box is periodic. Other boundaries are set with
``domain_boundaries``: see :doc:`/user_guide/boundary_conditions`.


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

  Fix: use the exact names above, and a combination from the table (or from
  the list printed after the message). The Python side does not check the
  names, so the same message appears for ``riemann="rusanov"``, or for a
  limiter with a non-Linear reconstruction such as
  ``reconstruction="WENOZ", limiter="VanLeer"``
  (``1,WENOZ,VanLeer,Rusanov,false``). Leaving ``limiter`` out gives the same
  message with an empty name: write ``limiter="None"`` for the schemes that
  take no limiter.

* **``nu`` without the Hall term**:

  .. code-block:: text

     RuntimeError: Error - hyper-resistivity (nu > 0) requires the Hall term to be enabled

  Fix: add ``hall=True`` (and a combination built with ``hall`` true, such as
  1D ``Linear``/``VanLeer``/``Rusanov``), or set ``nu=0``.

* **``resistivity`` or ``hyper_resistivity`` in an MHD run.** Both are hybrid
  options: an MHD run accepts them and ignores them. Use ``eta`` and ``nu``.

* **No time stepper**:

  .. code-block:: text

     ValueError: Error: MHDModel requires a non-empty mhd_timestepper

* **Unknown time stepper**, for example ``"rk4"``:

  .. code-block:: text

     ValueError: TimeIntegratorType: unknown value 'rk4', expected one of ['euler', 'tvdrk2', 'tvdrk3', 'ssprk4_5']

* **``gamma`` of 1 or less**:

  .. code-block:: text

     ValueError: 'gamma' must be a scalar greater than 1, got 1.0

* **``max_mhd_level`` too large**:

  .. code-block:: text

     ValueError: Error: max_mhd_level(3) should be less or equal to max_nbr_levels(1)

* **Using the hybrid model in an MHD run.** ``MaxwellianFluidModel`` belongs to the
  hybrid model. In a run declared with ``model_options=["MHDModel"]`` it fails
  with an unhelpful message when the simulation starts:

  .. code-block:: text

     ValueError: HybridGridLayoutFor: interp_order is 0 (no hybrid model) and no explicit ghosts_nbr was given

  Fix: use :py:class:`~pyphare.pharein.MHDModel`.

* **Reading results with the wrong accessor**: MHD quantities are read with
  ``run.GetMHDrho``, ``GetMHDV``, ``GetMHDP`` and the magnetic field with
  ``run.GetB`` (write it with ``ElectromagDiagnostics(quantity="B", ...)``);
  ``GetMHDrhoV`` and ``GetMHDEtot`` read ``rhoV`` and ``Etot``. All of them
  accept ``merged=True``, like ``GetB``.
* **Hybrid options in an MHD run.** ``interp_order`` is accepted by Python, but
  it asks for a hybrid+MHD module that is not built: the run stops with
  ``No module named 'pybindlibs.cpp_1_1_2_...'``. ``refined_particle_nbr``
  alone is rejected at once: ``Error: refined_particle_nbr is meaningless
  without a hybrid model (no interp_order given), it should be omitted``.
  Leave both out.


Going further
-------------

* :doc:`/user_guide/boundary_conditions`: open, reflective and inflow boundaries.
* :doc:`/user_guide/diagnostics`: choosing what to write and when.
* :doc:`/numerics/mhd_finite_volume`: how the scheme works.
* :doc:`/theory/mhd`: the equations.
