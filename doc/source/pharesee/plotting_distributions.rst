
===============================
Plotting particle distributions
===============================

PHARE stores particles as raw arrays of cell indices and fractional offsets
rather than as physical coordinates.  This page explains how to read those
arrays, compute physical positions and velocities, and produce the two most
common analysis plots: phase-space scatter plots and velocity-distribution
histograms.

.. contents:: On this page
   :local:
   :depth: 2


Loading particle data
---------------------

Use :meth:`~pyphare.pharesee.run.Run.GetParticles` to load particle data for
a named ion population at a given simulation time.  The return value is a
:class:`~pyphare.pharesee.hierarchy.PatchHierarchy` whose patch data hold
:class:`~pyphare.pharesee.particles.Particles` objects instead of field
arrays.

.. code-block:: python

   from pyphare.pharesee.run import Run

   run = Run("diags")
   particles = run.GetParticles(0.5, "protons")

To load several populations at once, pass a list or tuple of names.  All
populations are merged into a single hierarchy:

.. code-block:: python

   hier = run.GetParticles(0.5, ["protons", "alphas"])


Particle data structure
-----------------------

Each :class:`~pyphare.pharesee.particles.Particles` instance exposes the
following attributes.

.. list-table::
   :header-rows: 1
   :widths: 15 20 65

   * - Attribute
     - Shape
     - Description
   * - ``iCells``
     - ``[N, dim]``
     - Integer cell index of each particle (AMR cell coordinates at the
       patch level).
   * - ``deltas``
     - ``[N, dim]``
     - Fractional offset within the cell, in the half-open interval
       ``[0, 1)``.
   * - ``v``
     - ``[N, 3]``
     - Velocity vector.  Always three components even in 1-D and 2-D
       simulations; unused components are still present in the array.
   * - ``weights``
     - ``[N]``
     - Statistical weight of each macro-particle.  Use these as histogram
       weights to recover physical-unit distribution functions.
   * - ``charges``
     - ``[N]``
     - Charge of each macro-particle (in simulation units).
   * - ``dl``
     - ``[N, dim]``
     - Grid spacing at the particle's location.  Each row duplicates the
       cell width of the patch the particle belongs to.

Computing physical position
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The :class:`~pyphare.pharesee.particles.Particles` class provides ``x`` and
``y`` properties that compute the physical coordinate directly:

.. code-block:: python

   p = pd.dataset           # Particles object from a patch's patch_data
   x = p.x                  # dl[:, 0] * (iCells[:, 0] + deltas[:, 0])
   y = p.y                  # dl[:, 1] * (iCells[:, 1] + deltas[:, 1])  (2-D only)

For other axes, or when you need to be explicit, call ``p.xyz(i)``
directly:

.. code-block:: python

   import numpy as np
   x = p.dl[:, 0] * (p.iCells[:, 0] + p.deltas[:, 0])


Phase-space plots
-----------------

A phase-space scatter plot shows particle position against one velocity
component.  Iterating over patches collects particles from every AMR patch
on the requested level.

.. code-block:: python

   import numpy as np
   import matplotlib.pyplot as plt
   from pyphare.pharesee.run import Run

   run = Run("diags")
   particles = run.GetParticles(0.5, "protons")

   fig, ax = plt.subplots(figsize=(8, 5))

   for patch in particles.level(0).patches:
       pd = list(patch.patch_datas.values())[0]
       p  = pd.dataset
       ax.scatter(p.x, p.v[:, 0], s=0.5, alpha=0.2, c="steelblue",
                  linewidths=0)

   ax.set_xlabel("x")
   ax.set_ylabel("vx")
   ax.set_title("Phase space (protons, t = 0.5)")
   fig.tight_layout()
   plt.show()

.. note::

   ``patch.patch_datas`` is a dict whose key is the population name.
   ``list(...values())[0]`` fetches the first (and usually only) entry.
   If you loaded multiple populations, iterate the keys explicitly.


1-D velocity distributions
---------------------------

Weighted histograms recover the physical velocity distribution function.
Collect velocity and weight arrays from all patches first, then pass them to
:func:`numpy.histogram` or :func:`matplotlib.pyplot.hist`.

.. code-block:: python

   import numpy as np
   import matplotlib.pyplot as plt
   from pyphare.pharesee.run import Run

   run = Run("diags")
   particles = run.GetParticles(0.5, "protons")

   all_vx      = []
   all_weights = []

   for patch in particles.level(0).patches:
       p = list(patch.patch_datas.values())[0].dataset
       all_vx.append(p.v[:, 0])
       all_weights.append(p.weights)

   all_vx      = np.concatenate(all_vx)
   all_weights = np.concatenate(all_weights)

   fig, ax = plt.subplots()
   ax.hist(all_vx, bins=100, weights=all_weights, density=True,
           histtype="step", color="steelblue")
   ax.set_xlabel("vx")
   ax.set_ylabel("f(vx)")
   ax.set_title("Velocity distribution (protons, t = 0.5)")
   fig.tight_layout()
   plt.show()

.. tip::

   Always pass ``weights=all_weights`` to recover the physical distribution.
   An unweighted histogram treats all macro-particles as equal regardless of
   their physical content, which is generally wrong.


2-D velocity-space distributions
----------------------------------

A 2-D histogram in velocity space reveals temperature anisotropy, beams, and
ring distributions that are invisible in 1-D marginals.

.. code-block:: python

   import numpy as np
   import matplotlib.pyplot as plt
   from pyphare.pharesee.run import Run

   run = Run("diags")
   particles = run.GetParticles(0.5, "protons")

   all_vx, all_vy, all_weights = [], [], []

   for patch in particles.level(0).patches:
       p = list(patch.patch_datas.values())[0].dataset
       all_vx.append(p.v[:, 0])
       all_vy.append(p.v[:, 1])
       all_weights.append(p.weights)

   vx = np.concatenate(all_vx)
   vy = np.concatenate(all_vy)
   w  = np.concatenate(all_weights)

   fig, ax = plt.subplots(figsize=(6, 5))
   h, xedges, yedges, im = ax.hist2d(
       vx, vy, bins=80, weights=w, density=True, cmap="inferno"
   )
   fig.colorbar(im, ax=ax, label="f(vx, vy)")
   ax.set_xlabel("vx")
   ax.set_ylabel("vy")
   ax.set_aspect("equal")
   ax.set_title("Velocity-space distribution (protons, t = 0.5)")
   fig.tight_layout()
   plt.show()


Convenience function ``dist_plot``
------------------------------------

:func:`pyphare.pharesee.plotting.dist_plot` automates weighted 2-D histogram
plots for any combination of position and velocity axes.  It accepts a single
:class:`~pyphare.pharesee.particles.Particles` object, a list, or a dict of
populations (which it aggregates automatically via
:func:`~pyphare.pharesee.particles.aggregate`).

.. code-block:: python

   from pyphare.pharesee.run import Run
   from pyphare.pharesee.particles import aggregate
   from pyphare.pharesee.plotting import dist_plot

   run = Run("diags")
   hier = run.GetParticles(0.5, "protons")

   # collect Particles from all patches on level 0
   parts = aggregate(
       [list(p.patch_datas.values())[0].dataset
        for p in hier.level(0).patches]
   )

   fig, ax = dist_plot(parts, axis=("Vx", "Vy"), bins=(80, 80),
                        cmap="jet", color_scale="log",
                        xlabel="vx", ylabel="vy",
                        title="protons vx–vy at t = 0.5")
   fig.tight_layout()
   fig.savefig("dist_vx_vy.png", dpi=150)

Supported ``axis`` pairs:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - ``axis`` value
     - Description
   * - ``("Vx", "Vy")`` (default)
     - In-plane velocity space.
   * - ``("Vx", "Vz")``
     - Out-of-plane velocity components.
   * - ``("Vy", "Vz")``
     - Transverse velocity space.
   * - ``("x", "Vx")``
     - Phase space along *x* versus *vx*.
   * - ``("x", "Vy")``
     - Phase space along *x* versus *vy*.
   * - ``("x", "Vz")``
     - Phase space along *x* versus *vz*.

All keyword arguments accepted by ``dist_plot``:

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Keyword
     - Default
     - Description
   * - ``axis``
     - ``("Vx", "Vy")``
     - Pair of axes to plot (see table above).
   * - ``bins``
     - ``(50, 50)``
     - Number of histogram bins in each dimension.
   * - ``cmap``
     - ``"jet"``
     - Matplotlib colormap name.
   * - ``color_scale``
     - ``"log"``
     - ``"log"`` or ``"linear"`` normalisation for the colour axis.
   * - ``color_min`` / ``color_max``
     - histogram min / max
     - Override the colour-axis limits.
   * - ``gaussian_filter_sigma``
     - (not applied)
     - Sigma for :func:`scipy.ndimage.gaussian_filter` applied to the
       histogram image.  Cannot be combined with ``median_filter_size``.
   * - ``median_filter_size``
     - (not applied)
     - Size for :func:`scipy.ndimage.median_filter`.  Cannot be combined
       with ``gaussian_filter_sigma``.
   * - ``kde``
     - ``False``
     - Overlay kernel-density-estimate contours (requires *seaborn*).
   * - ``bulk``
     - ``False``
     - Draw dashed lines at the weighted bulk velocity along each axis.
   * - ``ax``
     - new axes
     - Existing :class:`matplotlib.axes.Axes` to draw into.
   * - ``xlabel`` / ``ylabel``
     - axis name
     - Axis labels.
   * - ``title``
     - ``""``
     - Axes title.
   * - ``xlim`` / ``ylim``
     - (not set)
     - Axis limits.
   * - ``filename``
     - (not saved)
     - If given, save the figure to this path.
   * - ``interp``
     - ``False``
     - If ``True``, also return a
       :class:`scipy.interpolate.LinearNDInterpolator` of the distribution,
       plus the bin-centre arrays.  Return value is then
       ``(fig, ax, interp_fn, xbins, ybins)``.

The standard return value is ``(fig, ax)``.


Spatial sub-selection before plotting
--------------------------------------

For large simulations the full particle set can be expensive to plot.
Restrict the spatial region either at load time (cheapest) or afterwards
with :meth:`~pyphare.pharesee.particles.Particles.select`.

Load-time restriction via ``selection_box``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Load only particles in the sub-domain x in [20, 40] for a 1-D run
   from pyphare.core.box import Box

   hier = run.GetParticles(0.5, "protons",
                           selection_box=Box([20], [40]))

Post-load restriction with ``Particles.select``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from pyphare.core.box import Box

   p_all = list(hier.level(0).patches[0].patch_datas.values())[0].dataset
   p_sub = p_all.select(Box([20], [40]), box_type="pos")

``box_type="pos"`` selects by physical position; ``box_type="cell"`` selects
by AMR cell index (the default).


Tips
----

- **Always use weights** when computing histograms.  Skipping them produces
  an unphysical result whenever macro-particle weights are not uniform.
- **Downsample for scatter plots** on large datasets.  A scatter plot of
  10 million particles is slow and visually indistinguishable from 100 000:

  .. code-block:: python

     rng  = np.random.default_rng(42)
     mask = rng.integers(0, len(vx), size=100_000)
     ax.scatter(x[mask], vx[mask], s=0.5, alpha=0.1)

- **Aggregate across patches** before passing to ``dist_plot``.  The
  function accepts a list of :class:`~pyphare.pharesee.particles.Particles`
  and calls :func:`~pyphare.pharesee.particles.aggregate` internally, so you
  can also pass the list directly.
- **Multi-level data** requires care.  Finer AMR levels overlap coarser ones;
  plotting all levels together counts particles twice.  Either restrict to a
  single level or apply a spatial mask to exclude the coarse-level overlap
  regions before aggregating.
