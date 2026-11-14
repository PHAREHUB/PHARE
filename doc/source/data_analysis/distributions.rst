===============================
Plotting particle distributions
===============================

**Goal:** write the macro-particles of a run to disk, read them back, and plot
a velocity distribution or a phase-space map (position against velocity).

.. contents:: On this page
   :local:
   :depth: 1


Minimal example
---------------

Take ``first_simulation.py`` from :doc:`/getting_started/first_simulation` and
add one diagnostic block before ``Simulator(sim).run()``:

.. code-block:: python

   times = [0.0, 25.0, 50.0]
   ph.ParticleDiagnostics(quantity="domain", population_name="protons",
                          write_timestamps=times)

``population_name`` is the name of the population in your
``MaxwellianFluidModel`` (here ``protons``). After the run, read the
particles at t = 25 and plot the distribution of the velocity components
:math:`v_x` and :math:`v_y`:

.. code-block:: python

   from pyphare.pharesee.run import Run

   run = Run("first_run")
   particles = run.GetParticles(25.0, "protons")      # a patch hierarchy

   final, (fig, ax) = particles.dist_plot(
       time=25.0, axis=("Vx", "Vy"), bins=(60, 60),
       title="protons at t = 25", filename="vxvy.png",
   )

The file ``first_run/ions_pop_protons_domain.h5`` holds the particles.
``GetParticles`` returns a patch hierarchy (see :doc:`/data_analysis/patch_hierarchy`)
whose data are lists of particles, instead of fields.

.. important:: Particles are written only with the ``"phareh5"`` format. With
   ``"pharevtkhdf"``, a ``ParticleDiagnostics`` block is accepted but **nothing
   is written** (no file, no message). Also expect large files: every particle
   is stored.


What it does
------------

``dist_plot`` collects the particles of the chosen levels, makes a histogram
and draws it. It returns two things: a dictionary of the particles used
(``{"protons_domain": Particles}``) and ``(fig, ax)`` of the figure. A
``Particles`` object has ``x`` (positions), ``v`` (an array of shape
``(N, 3)``: the three velocity components), ``weights`` and ``size()`` (number
of particles).

.. code-block:: python

   p = final["protons_domain"]
   print(p.size(), p.v.shape)       # 40000 (40000, 3)
   print(p.v[:, 0].mean())          # mean of vx

Each particle counts with its weight in the histogram.

At t = 0, the initial condition gives :math:`v_y = -0.1\cos(kx)` plus a
thermal spread, so a map of :math:`v_y` against :math:`x` shows a band that
follows the cosine (at later times the band has moved with the wave). The figure
below is drawn from the formula, not from a PHARE run: your map has the same
shape, with the particle noise on top.

.. plot::
   :caption: Expected shape of ``axis=("x", "Vy")`` for the first simulation at
             t = 0 (formula and random numbers, not a PHARE run): the bulk
             velocity follows :math:`-0.1\cos(kx)`, with a thermal spread of 0.1.

   import numpy as np
   import matplotlib.pyplot as plt

   plt.rcParams["font.size"] = 11
   rng = np.random.default_rng(0)
   L, a, vth = 100.0, 0.1, 0.1
   x = rng.uniform(0, L, 40000)
   vy = -a * np.cos(2 * np.pi * x / L) + vth * rng.normal(size=x.size)
   fig, ax = plt.subplots(figsize=(6, 3.4))
   ax.hist2d(x, vy, bins=(100, 60), cmap="jet", cmin=1)
   ax.set_xlabel("x")
   ax.set_ylabel("Vy")
   fig.tight_layout()

Other maps from the same data:

.. code-block:: python

   # position against Vy
   particles.dist_plot(time=25.0, axis=("x", "Vy"), bins=(100, 60),
                       filename="xvy.png")

   # one velocity component: a curve
   particles.dist_plot(time=25.0, axis="Vx", bins=60, filename="vx.png")

   # smooth the map, linear colour scale
   particles.dist_plot(time=25.0, axis=("Vy", "Vz"), color_scale="linear",
                       cmap="viridis", gaussian_filter_sigma=(1, 1),
                       filename="vyvz.png")

To keep only some particles, pass a function to ``select`` that takes a
``Particles`` object and returns the selected ones; here the particles with
:math:`0 \le x \le 50`:

.. code-block:: python

   from pyphare.core.box import Box

   final, _ = particles.dist_plot(
       time=25.0, axis=("Vx", "Vy"),
       select=lambda p: p.select(Box(0, 50), box_type="pos"),
   )
   print(final["protons_domain"].v.shape)    # about half the particles


Options
-------

Options of ``hierarchy.dist_plot(**kwargs)`` (the ones that choose the particles):

.. list-table::
   :header-rows: 1
   :widths: 22 56 22

   * - Option
     - Meaning
     - Default
   * - ``time``
     - time to plot
     - first time of the hierarchy
   * - ``levels``
     - refinement levels whose particles are used
     - ``(0,)``
   * - ``finest``
     - ``True``: use, at each place, the particles of the finest level
       available (no double counting); ``levels`` is then ignored
     - ``False``
   * - ``pop``
     - names of populations to use (list); the hierarchy keys are the file
       names, for example ``"protons_domain"``
     - all populations in the hierarchy
   * - ``select``
     - function ``Particles -> Particles`` applied before plotting
     - none
   * - ``vmin``, ``vmax``, ``dv``
     - velocity bin range and width, used only with ``finest=True``
     - -2, 2, 0.05

Options of the drawing itself (passed on to the plotting function):

.. list-table::
   :header-rows: 1
   :widths: 22 56 22

   * - Option
     - Meaning
     - Default
   * - ``axis``
     - ``("Vx", "Vy")`` and the other pairs of ``Vx``, ``Vy``, ``Vz``;
       ``("x", "Vx")``, ``("x", "Vy")``, ``("x", "Vz")``; or one velocity name for a
       curve
     - ``("Vx", "Vy")``
   * - ``bins``
     - number of bins, one number per axis (one number for a curve)
     - ``(50, 50)``, 50
   * - ``cmap``
     - colour map
     - ``"jet"``
   * - ``color_scale``
     - ``"log"`` or ``"linear"``
     - ``"log"``
   * - ``color_min``, ``color_max``
     - colour range (``color_min`` is at least 1e-4)
     - data min and max
   * - ``gaussian_filter_sigma`` / ``median_filter_size``
     - smooth the map; the two cannot be used together
     - none
   * - ``cuts``
     - for ``axis="Vx"`` only: position range ``[(xmin, xmax)]`` of the
       particles to keep
     - none
   * - ``title``, ``xlabel``, ``ylabel``, ``xlim``, ``ylim``
     - text and axis limits
     - axis names
   * - ``ax``
     - matplotlib axes to draw on
     - a new figure
   * - ``filename``
     - save the figure under this name
     - not saved
   * - ``kde``
     - ``True`` adds contours of a density estimate (needs ``seaborn``)
     - ``False``

``ParticleDiagnostics`` has two options: ``quantity`` (``"domain"``: the
particles inside the patches, the right choice for analysis; ``"levelGhost"``:
particles in the ghost region between levels, for debugging refinement) and
``population_name``. The list of all diagnostics is in
:doc:`/user_guide/diagnostics`.

To plot only some populations, pass the names the hierarchy uses, for example
``pop=["protons_domain"]``. For a run with several populations, read them all
at once with ``run.GetParticles(t, ["pop1", "pop2"])``.


Common mistakes
---------------

* ``FileNotFoundError: ... 'first_run/ions_pop_protons_domain.h5' ...`` — no
  particle diagnostic was requested for that population, or it was written
  with ``"pharevtkhdf"`` (nothing is written there). Use ``"phareh5"``.
* ``KeyError: 'protons'`` when reading ``final`` — the dictionary key is
  ``"protons_domain"`` (population and diagnostic), not the population name.
* An ``AssertionError`` with no message — you called ``dist_plot`` on a field
  hierarchy (``GetB``...) instead of the output of ``GetParticles``.

Next: :doc:`plotting_fields` for fields, :doc:`reading_data` for the other getters.
