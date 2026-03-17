
================
Plotting fields
================

PHARE diagnostic output stores field data (magnetic field, electric field,
current density, etc.) as HDF5 files.  The ``pharesee`` analysis layer reads
those files back into a :class:`~pyphare.pharesee.hierarchy.PatchHierarchy`
that mirrors the AMR patch structure used during the simulation.  This page
shows how to turn that hierarchy into plots with matplotlib.

.. contents:: On this page
   :local:
   :depth: 2


1D field plotting
-----------------

For a 1D simulation, load the run and request the magnetic field at a
specific simulation time.  Passing ``merged=True`` interpolates all patches
onto a single uniform grid so that each AMR level looks like one contiguous
array — the easiest way to get a quick plot.

.. code-block:: python

   from pyphare.pharesee.run import Run
   import matplotlib.pyplot as plt

   run = Run("diags")

   # merged=True merges all patches of each level onto a uniform grid.
   # interp="nearest" is the default interpolation mode.
   B = run.GetB(0.5, merged=True, interp="nearest")

   # For 1D merged data each level is represented by a single patch.
   for patch in B.level(0).patches:
       by = patch.patch_datas["By"]
       plt.plot(by.x, by.dataset)

   plt.xlabel("x")
   plt.ylabel("By")
   plt.title("By at t = 0.5")
   plt.show()

.. note::

   ``run.GetB()`` returns a :class:`~pyphare.pharesee.run.VectorField` whose
   ``Bx``, ``By``, and ``Bz`` sub-hierarchies are accessible as attributes.
   You can pass the component hierarchy directly to plotting routines:
   ``B.By.plot1d(qty="By")``.

When ``merged=False`` (the default) the hierarchy contains one patch per
original AMR patch.  Ghost cells are present but excluded automatically by
``plot1d`` when it accesses data via ``pd[patch.box]``.


Plotting time evolution
-----------------------

Plot the same component at several times on the same axes to show how a field
evolves.

.. code-block:: python

   import matplotlib.pyplot as plt
   from pyphare.pharesee.run import Run

   run = Run("diags")
   fig, ax = plt.subplots()

   for t in [0.0, 0.25, 0.5, 0.75, 1.0]:
       B = run.GetB(t, merged=True)
       for patch in B.By.level(0).patches:
           by = patch.patch_datas["By"]
           ax.plot(by.x, by.dataset, label=f"t = {t:.2f}")

   ax.set_xlabel("x")
   ax.set_ylabel("By")
   ax.legend()
   fig.tight_layout()
   plt.show()


2D field plotting
-----------------

For 2D simulations use :meth:`~pyphare.pharesee.hierarchy.PatchHierarchy.plot2d`
or construct a ``pcolormesh`` manually.  Both approaches handle the fact that
there may be multiple AMR patches covering the domain.

Using the convenience method
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from pyphare.pharesee.run import Run
   import matplotlib.pyplot as plt

   run = Run("diags")
   Bz = run.GetB(0.5).Bz   # single-component hierarchy

   fig, ax = plt.subplots(figsize=(8, 6))
   Bz.plot2d(qty="Bz", ax=ax, cmap="RdBu_r", title="Bz at t = 0.5")
   plt.show()

Accepted keyword arguments for ``plot2d``:

.. list-table::
   :header-rows: 1
   :widths: 20 20 60

   * - Keyword
     - Default
     - Description
   * - ``qty``
     - auto (if unique)
     - Name of the quantity to plot.
   * - ``time``
     - first available time
     - Simulation time to visualise.
   * - ``levels``
     - all levels
     - List of AMR level numbers to include.
   * - ``cmap``
     - ``"Spectral_r"``
     - Matplotlib colormap name.
   * - ``vmin`` / ``vmax``
     - global min/max
     - Color scale limits.
   * - ``plot_patches``
     - ``False``
     - Overlay AMR patch boundaries (see below).
   * - ``patchcolors``
     - ``"k"`` for all levels
     - Dict or list of colors per level for patch outlines.
   * - ``lw``
     - ``1`` for all levels
     - Line width(s) for patch outlines.
   * - ``ls``
     - ``"-"`` for all levels
     - Line style(s) for patch outlines.
   * - ``aspect``
     - ``"equal"``
     - Axes aspect ratio.
   * - ``xlabel`` / ``ylabel``
     - ``"x"`` / ``"y"``
     - Axis labels.
   * - ``xlim`` / ``ylim``
     - (not set)
     - Axis limits.
   * - ``ax``
     - new axes
     - Existing :class:`matplotlib.axes.Axes` to draw into.
   * - ``filename``
     - (not saved)
     - If given, save the figure to this path.
   * - ``dpi``
     - ``200``
     - Resolution when saving.

Manual ``pcolormesh``
~~~~~~~~~~~~~~~~~~~~~

When you need full control over the colour mapping or want to combine data
from several quantities on the same axes, iterate patches directly:

.. code-block:: python

   import numpy as np
   import matplotlib.pyplot as plt
   from pyphare.pharesee.run import Run

   run = Run("diags")
   Bz = run.GetB(0.5).Bz

   fig, ax = plt.subplots(figsize=(8, 6))

   for patch in Bz.level(0).patches:
       pdat = patch.patch_datas["Bz"]
       ng = pdat.ghosts_nbr
       data = pdat.dataset[ng[0]:-ng[0], ng[1]:-ng[1]]
       x = pdat.x[ng[0]:-ng[0]]
       y = pdat.y[ng[1]:-ng[1]]
       dx, dy = pdat.layout.dl
       # pcolormesh needs cell *edges*, not cell centres
       x_edges = np.append(x - dx * 0.5, x[-1] + dx * 0.5)
       y_edges = np.append(y - dy * 0.5, y[-1] + dy * 0.5)
       ax.pcolormesh(x_edges, y_edges, data.T, cmap="RdBu_r")

   ax.set_aspect("equal")
   ax.set_xlabel("x")
   ax.set_ylabel("y")
   plt.colorbar(ax.collections[0], ax=ax)
   plt.show()


Contour plots and magnetic flux
--------------------------------

The vector potential :math:`A_z` is available via
:meth:`~pyphare.pharesee.run.Run.GetMagneticFlux`.  Its contours correspond
to magnetic field lines, making it the standard overlay for 2D reconnection
diagnostics.

.. code-block:: python

   import numpy as np
   import matplotlib.pyplot as plt
   from pyphare.pharesee.run import Run

   run = Run("diags")

   # GetMagneticFlux returns (Az_array, (xn, yn)) where xn, yn are the
   # uniform node-centred grids used for interpolation.
   Az, (xn, yn) = run.GetMagneticFlux(0.5)

   Bz = run.GetB(0.5).Bz

   fig, ax = plt.subplots(figsize=(8, 6))
   Bz.plot2d(qty="Bz", ax=ax, cmap="RdBu_r")
   ax.contour(xn, yn, Az.T, levels=20, colors="k", linewidths=0.6)
   ax.set_title("Bz with magnetic field lines (t = 0.5)")
   plt.show()


Convenience methods ``plot1d`` and ``plot2d``
---------------------------------------------

:class:`~pyphare.pharesee.hierarchy.PatchHierarchy` exposes both methods
directly.  The dimension-agnostic entry point ``hier.plot()`` dispatches to
the correct one automatically.

``plot1d(**kwargs)``
~~~~~~~~~~~~~~~~~~~~~

Plots every patch on the requested levels as separate lines, labelled
``L{level}P{patch}``.

.. code-block:: python

   B = run.GetB(0.5).By
   B.plot1d(qty="By", levels=(0, 1), ls="-", color="steelblue",
            xlabel="x", ylabel="By", title="By — levels 0 and 1")

Accepted keyword arguments:

.. list-table::
   :header-rows: 1
   :widths: 20 20 60

   * - Keyword
     - Default
     - Description
   * - ``qty``
     - auto (if unique)
     - Quantity name to plot.
   * - ``time``
     - first available time
     - Simulation time to visualise.
   * - ``levels``
     - ``(0,)``
     - Tuple of AMR level numbers to include.
   * - ``ax``
     - new axes
     - Existing :class:`matplotlib.axes.Axes` to draw into.
   * - ``marker``
     - ``""``
     - Matplotlib marker string.
   * - ``ls``
     - ``"--"``
     - Line style.
   * - ``color``
     - ``"k"``
     - Line colour.
   * - ``xlabel`` / ``ylabel``
     - ``"x"`` / quantity name
     - Axis labels.
   * - ``title``
     - ``""``
     - Axes title.
   * - ``xlim`` / ``ylim``
     - (not set)
     - Axis limits.
   * - ``legend``
     - (not shown)
     - Pass any non-``None`` value to enable the legend.
   * - ``filename``
     - (not saved)
     - If given, save the figure to this path.

``plot2d(**kwargs)``
~~~~~~~~~~~~~~~~~~~~~

See the keyword table in the *2D field plotting* section above.  The
``plot2d`` method also returns ``(fig, ax)`` so you can continue decorating
the figure.

.. code-block:: python

   fig, ax = B.Bz.plot2d(qty="Bz", plot_patches=True,
                          patchcolors={0: "k", 1: "r"},
                          title="Bz with AMR patch boundaries")


AMR patch boundaries
--------------------

Overlaying patch boundaries helps verify that the refinement regions are
placed where expected.

Using the built-in ``plot_patches=True`` flag
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Pass ``plot_patches=True`` to ``plot2d`` (see the keyword table above).  Each
refined level gets its own colour:

.. code-block:: python

   fig, ax = Bz.plot2d(
       qty="Bz",
       plot_patches=True,
       patchcolors={0: "grey", 1: "red", 2: "orange"},
       lw={0: 1, 1: 2, 2: 2},
   )

Drawing rectangles manually
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For 2D simulations the patch ``box`` is in *cell* index space.  Multiply by
the cell size ``dl`` to get physical coordinates:

.. code-block:: python

   import matplotlib.pyplot as plt
   import matplotlib.patches as mpatches
   from pyphare.pharesee.run import Run

   run = Run("diags")
   Bz = run.GetB(0.5).Bz

   fig, ax = plt.subplots(figsize=(8, 6))
   Bz.plot2d(qty="Bz", ax=ax, cmap="RdBu_r")

   colors = {0: "black", 1: "red"}
   for ilvl, lvl in Bz.levels().items():
       for patch in lvl.patches:
           dl = patch.dl
           lo = patch.box.lower * dl
           shape = patch.box.shape * dl
           rect = mpatches.Rectangle(
               lo, shape[0], shape[1],
               linewidth=1.5,
               edgecolor=colors.get(ilvl, "blue"),
               facecolor="none",
               linestyle="--",
           )
           ax.add_patch(rect)

   plt.show()


Multi-panel figures
-------------------

Plot all three components of a vector field side by side using
:func:`matplotlib.pyplot.subplots`.

.. code-block:: python

   import matplotlib.pyplot as plt
   from pyphare.pharesee.run import Run

   run = Run("diags")
   B = run.GetB(0.5)

   components = [("Bx", B.Bx), ("By", B.By), ("Bz", B.Bz)]
   fig, axes = plt.subplots(1, 3, figsize=(15, 5))

   for ax, (name, hier) in zip(axes, components):
       hier.plot2d(qty=name, ax=ax, title=name)

   fig.suptitle("Magnetic field components at t = 0.5")
   fig.tight_layout()
   plt.show()

For a 1D simulation, replace ``plot2d`` with ``plot1d``:

.. code-block:: python

   fig, axes = plt.subplots(3, 1, figsize=(8, 9), sharex=True)

   for ax, (name, hier) in zip(axes, components):
       hier.plot1d(qty=name, ax=ax, ylabel=name)

   axes[-1].set_xlabel("x")
   fig.tight_layout()
   plt.show()
