===============
Plotting fields
===============

**Goal:** plot a field (magnetic field, density...) from a PHARE run, in 1D and
in 2D, choose which refinement levels to draw, and save the figure.

The 1D examples use the run of :doc:`reading_data` (the first simulation with
more diagnostics). The 2D examples use the small run described in
`A small 2D run`_.

.. contents:: On this page
   :local:
   :depth: 1


Minimal example (1D)
--------------------

.. code-block:: python

   import matplotlib.pyplot as plt
   from pyphare.pharesee.run import Run

   run = Run("first_run")

   # the magnetic field: merged=True gives a function you can plot yourself
   fig, ax = plt.subplots()
   for t in run.times("B"):
       by, (x,) = run.GetB(t, merged=True)["By"]
       ax.plot(x, by(x), label=f"t = {t:g}")
   ax.set_xlabel("x")
   ax.set_ylabel("By")
   ax.legend()
   fig.savefig("by.png")

   # the density: use the plot method of the patch hierarchy
   Ni = run.GetNi(25.0)
   Ni.plot(qty="value", time=25.0, filename="ni.png")

Two ways to plot, depending on what the getter returns (see
:doc:`reading_data`):

* ``merged=True`` gives, for each component, a function and its grid points.
  You draw with ``matplotlib`` as usual. It currently works only for ``GetB`` and
  ``GetE``.
* Without it, you get a patch hierarchy (:doc:`/data_analysis/patch_hierarchy`)
  that has its own ``plot`` method. It works for every quantity.

.. note:: **Known issue.** ``merged=True`` currently works only for ``GetB`` and
   ``GetE``. For the others (``GetNi``, ``GetVi``...) use the hierarchy
   ``plot`` method.


The ``plot`` method
-------------------

``hierarchy.plot(**kwargs)`` calls ``plot1d`` for a 1D run and ``plot2d`` for a
2D run (nothing is drawn for 3D). You always pass options as keywords.

Which ``qty``? A hierarchy can hold several quantities: the components of a
vector field are called ``"x"``, ``"y"``, ``"z"`` (so ``B.plot(qty="y")`` draws
By), and a scalar field has a single quantity called ``"value"``.

.. code-block:: python

   B = run.GetB(25.0)
   B.plot(qty="y", time=25.0, title="By at t = 25", ls="-", color="C0")

To put several plots in one figure, create the axes yourself and pass ``ax``:

.. code-block:: python

   fig, ax = plt.subplots(2, 1, sharex=True)
   run.GetB(50.0).plot(qty="y", time=50.0, ax=ax[0], ylabel="By")
   run.GetNi(50.0).plot(qty="value", time=50.0, ax=ax[1], ylabel="Ni")
   fig.savefig("two_panels.png")

1D options (``plot1d``)
^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 22 56 22

   * - Option
     - Meaning
     - Default
   * - ``qty``
     - quantity to draw (see above); required if the hierarchy holds several
     - the only one
   * - ``time``
     - time to draw
     - first time of the hierarchy
   * - ``levels``
     - tuple of refinement levels to draw, for example ``(0, 1)``
     - ``(0,)``
   * - ``ax``
     - matplotlib axes to draw on
     - a new figure
   * - ``marker``, ``ls``, ``color``
     - matplotlib line style, applied to all patches
     - ``""``, ``"--"``, ``"k"``
   * - ``title``, ``xlabel``, ``ylabel``
     - text
     - ``""``, ``"x"``, the quantity
   * - ``xlim``, ``ylim``
     - axis limits, a pair
     - automatic
   * - ``legend``
     - any value other than ``None`` draws a legend (one label per patch,
       ``L0P0`` = level 0, patch 0)
     - no legend
   * - ``filename``
     - save the figure under this name
     - not saved

``plot1d`` draws the patches of the chosen levels one by one: where the levels
overlap, you see both. It returns nothing; keep your own ``fig, ax`` (with the
``ax`` option) if you need them afterwards.


2D plots
--------

A small 2D run
^^^^^^^^^^^^^^

To test the 2D calls, this input script runs a uniform plasma with a small
density perturbation and a small wave in B (40 by 40 cells, 5 time units,
about 10 seconds on a laptop). The setup is the same as in
:doc:`/getting_started/first_simulation`.

.. code-block:: python

   import numpy as np
   import pyphare.pharein as ph
   from pyphare.simulator.simulator import Simulator

   sim = ph.Simulation(
       cells=(40, 40),
       dl=(0.5, 0.5),
       time_step=0.01,
       final_time=5.0,
       diag_options={"format": "phareh5",
                     "options": {"dir": "run2d", "mode": "overwrite"}},
   )

   Lx = Ly = 20.0
   kx, ky = 2 * np.pi / Lx, 2 * np.pi / Ly

   def density(x, y): return 1.0 + 0.1 * np.cos(kx * x) * np.cos(ky * y)
   def bx(x, y):      return 1.0
   def by(x, y):      return 0.05 * np.sin(kx * x)
   def bz(x, y):      return 0.0
   def v0(x, y):      return 0.0
   def vth(x, y):     return 0.1

   ph.MaxwellianFluidModel(
       bx=bx, by=by, bz=bz,
       protons={"charge": 1, "density": density,
                "vbulkx": v0, "vbulky": v0, "vbulkz": v0,
                "vthx": vth, "vthy": vth, "vthz": vth,
                "nbr_part_per_cell": 100},
   )
   ph.ElectronModel(closure="isothermal", Te=0.1)

   times = [0.0, 5.0]
   ph.ElectromagDiagnostics(quantity="B", write_timestamps=times)
   ph.ElectromagDiagnostics(quantity="E", write_timestamps=times)
   ph.FluidDiagnostics(quantity="charge_density", write_timestamps=times)
   ph.FluidDiagnostics(quantity="bulkVelocity", write_timestamps=times)

   Simulator(sim).run()

(The magnetic field is chosen so that ``By`` depends on ``x`` only: that keeps
:math:`\nabla\cdot\mathbf{B}=0` at the start.)

Drawing a map
^^^^^^^^^^^^^

.. code-block:: python

   from pyphare.pharesee.run import Run

   run = Run("run2d")
   Ni = run.GetNi(5.0)

   # simplest call: a colour map of the density, saved as a file
   Ni.plot(qty="value", time=5.0, filename="ni2d.png", title="Ni at t = 5")

   # more control: choose colours and range, draw the patch edges
   fig, ax = Ni.plot(qty="value", time=5.0, cmap="viridis", vmin=0.8, vmax=1.2,
                     plot_patches=True)
   fig.savefig("ni2d_b.png", dpi=100)

   # three components of B side by side, with your own axes
   import matplotlib.pyplot as plt
   B = run.GetB(5.0)
   fig, axes = plt.subplots(3, 1, figsize=(5, 12))
   for ax, c in zip(axes, "xyz"):
       B.plot(qty=c, time=5.0, ax=ax, title="B" + c)
   fig.tight_layout()
   fig.savefig("b2d.png")

``plot`` returns ``(fig, ax)`` in 2D. A colour bar is always added to the right.

2D options (``plot2d``)
^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 22 56 22

   * - Option
     - Meaning
     - Default
   * - ``qty``
     - quantity to draw
     - the only one, if the hierarchy holds one
   * - ``time``
     - time to draw
     - first time of the hierarchy
   * - ``levels``
     - refinement levels to draw
     - all levels
   * - ``ax``
     - matplotlib axes to draw on
     - a new figure
   * - ``cmap``
     - matplotlib colour map
     - ``"Spectral_r"``
   * - ``vmin``, ``vmax``
     - colour range
     - minimum and maximum of the quantity over the whole hierarchy (widened by 1e-6)
   * - ``plot_patches``
     - ``True`` draws the edge of every patch
     - ``False``
   * - ``patchcolors``, ``lw``, ``ls``
     - edge colour, width and style of patches, one per level: a dict
       ``{level: value}`` or a list, **not** a single value (``lw=2`` stops
       with ``TypeError: 'int' object is not iterable``)
     - ``"k"``, 1, ``"-"``
   * - ``aspect``
     - matplotlib aspect of the axes
     - ``"equal"``
   * - ``title``, ``xlabel``, ``ylabel``
     - text
     - ``""``, ``"x"``, ``"y"``
   * - ``xlim``, ``ylim``
     - axis limits, a pair
     - automatic
   * - ``legend``
     - any value other than ``None`` calls ``ax.legend()``
     - none
   * - ``filename``, ``dpi``
     - save under this name, at this resolution
     - not saved, 200

Zoom with ``xlim`` and ``ylim``, for example
``Ni.plot(qty="value", time=5.0, xlim=(0, 10), ylim=(0, 10))``.

A map of a merged field (``GetB(t, merged=True)``) is made with
``numpy.meshgrid`` and ``pcolormesh``:

.. code-block:: python

   import numpy as np
   import matplotlib.pyplot as plt

   by, (x, y) = run.GetB(5.0, merged=True)["By"]
   X, Y = np.meshgrid(x, y, indexing="ij")   # the function wants 2D arrays
   fig, ax = plt.subplots()
   im = ax.pcolormesh(X, Y, by(X, Y), shading="auto")
   fig.colorbar(im)
   fig.savefig("by_merged.png")


Choosing levels
---------------

A hierarchy can have several refinement levels (finer grids where the
simulation needs them). Use ``levels`` to say which to draw. In 1D the default
is level 0 only, in 2D all levels: pass ``levels=(0,)`` to draw only the
coarsest grid in 2D. The 2D runs of this page have a single level (0), so
``levels=(0,)`` is the only valid choice for them.


Saving a figure
---------------

Two ways: pass ``filename="name.png"`` to ``plot`` (2D also accepts ``dpi``),
or keep the figure and call ``fig.savefig("name.png")`` yourself, as in the
examples above. On a machine without a display (a cluster node), put this
before importing ``pyplot``:

.. code-block:: python

   import matplotlib
   matplotlib.use("Agg")


Common mistakes
---------------

* ``ValueError: multiple quantities in patch, please specify a quantity in x y z``
  (1D) — you forgot ``qty`` on a vector field. Add ``qty="y"``.
* ``KeyError: None`` (2D, same cause) — add ``qty``.
* ``KeyError: 'By'`` — on a hierarchy the components are ``"x"``, ``"y"``,
  ``"z"``, not ``"Bx"``... (the names ``"Bx"``... are for ``merged=True``).
* ``KeyError: "Unable to synchronously open object (object '7.0000000000'
  doesn't exist)"`` — you asked for a time that was not written;
  check ``run.times("B")``.
* ``AttributeError: 'dict' object has no attribute 'domain_box'`` — you used
  ``merged=True`` with an ion quantity (known issue above).

Next: :doc:`distributions` for particles. Reading details are in :doc:`reading_data`.
