============================
Looking at your first result
============================

**Goal:** write more quantities from :doc:`your first simulation
<first_simulation>`, find them on disk, read them in Python, and check one more
piece of physics: in an Alfvén wave, the ion velocity follows the magnetic field.

.. contents:: On this page
   :local:
   :depth: 1


Step 1 — write more quantities
------------------------------

In ``first_simulation.py``, replace the line with ``ph.ElectromagDiagnostics``
by these five lines, then run the script again (about 15 seconds):

.. code-block:: python

   times = [0.0, 25.0, 50.0]
   ph.ElectromagDiagnostics(quantity="B", write_timestamps=times)
   ph.ElectromagDiagnostics(quantity="E", write_timestamps=times)
   ph.FluidDiagnostics(quantity="charge_density", write_timestamps=times)
   ph.FluidDiagnostics(quantity="bulkVelocity", write_timestamps=times)

``FluidDiagnostics`` writes *moments* of the ions: quantities computed from
all the particles in each cell, such as their density or their mean (bulk)
velocity. The full list of what PHARE can write is in
:doc:`/user_guide/diagnostics`.


Step 2 — what is on disk
------------------------

The folder ``first_run`` now holds one file per quantity:

.. list-table::
   :header-rows: 1
   :widths: 35 35 30

   * - File
     - Contains
     - Read with
   * - ``EM_B.h5``
     - magnetic field :math:`\mathbf{B}`
     - ``run.GetB(t)``
   * - ``EM_E.h5``
     - electric field :math:`\mathbf{E}`
     - ``run.GetE(t)``
   * - ``ions_charge_density.h5``
     - ion charge density (here equal to the number density: protons have
       charge 1)
     - ``run.GetNi(t)``
   * - ``ions_bulkVelocity.h5``
     - ion bulk velocity :math:`\mathbf{V}_i`
     - ``run.GetVi(t)``

Each file holds all the output times of its quantity. To see which times are
available:

.. code-block:: python

   from pyphare.pharesee.run import Run

   run = Run("first_run")
   print(run.all_times())
   # {'E': array([ 0., 25., 50.]), 'B': array([ 0., 25., 50.]),
   #  'Vi': array([ 0., 25., 50.]), 'Ni': array([ 0., 25., 50.])}


Step 3 — read and plot
----------------------

PHARE does not store a field as one big array. The domain is cut into
rectangular pieces of grid called **patches**; with mesh refinement, finer
patches sit on top of coarser ones, in *levels*. A getter such as
``run.GetB(t)`` returns this collection of patches (more in
:doc:`/data_analysis/patch_hierarchy`). The picture below shows an example
domain cut into four patches; your first run has a single one.

.. raw:: html

   <svg viewBox="0 0 600 120" style="max-width:600px;width:100%;height:auto" role="img"
        aria-label="A 1D domain from 0 to 100 cut into four patches of 25, each with ghost points just outside its edges">
     <style>
       .p { fill: #4c8bf5; fill-opacity: .18; stroke: #4c8bf5; stroke-width: 2 }
       .g { fill: #e8710a; }
       .t { font: 14px sans-serif; fill: currentColor }
     </style>
     <rect class="p" x="20"  y="30" width="140" height="40"/>
     <rect class="p" x="160" y="30" width="140" height="40"/>
     <rect class="p" x="300" y="30" width="140" height="40"/>
     <rect class="p" x="440" y="30" width="140" height="40"/>
     <text class="t" x="70"  y="55">patch 0</text>
     <text class="t" x="210" y="55">patch 1</text>
     <text class="t" x="350" y="55">patch 2</text>
     <text class="t" x="490" y="55">patch 3</text>
     <circle class="g" cx="10" cy="85" r="5"/><circle class="g" cx="590" cy="85" r="5"/>
     <circle class="g" cx="150" cy="85" r="5"/><circle class="g" cx="170" cy="85" r="5"/>
     <circle class="g" cx="290" cy="85" r="5"/><circle class="g" cx="310" cy="85" r="5"/>
     <circle class="g" cx="430" cy="85" r="5"/><circle class="g" cx="450" cy="85" r="5"/>
     <text class="t" x="20"  y="110">x = 0</text>
     <text class="t" x="520" y="110">x = 100</text>
     <text class="t" x="215" y="20">orange: ghost points, a copy of the neighbour's values</text>
   </svg>

Each patch also stores a few extra points just outside its edges, the *ghost*
points (orange), which it needs to compute derivatives. When you extract the
data, you drop them and keep the inside of each patch. This small function does
that for the coarsest level, for any quantity:

.. code-block:: python

   import numpy as np
   import matplotlib.pyplot as plt
   from pyphare.pharesee.run import Run


   def profile(field, component):
       """x and values of one component on the coarsest level, patch after patch."""
       xs, values = [], []
       for patch in field.level(0).patches:
           data = patch.patch_datas[component]
           n = data.ghosts_nbr[0]              # number of ghost points on each side
           xs.append(data.x[n:-n])             # positions, without ghosts
           values.append(data[patch.box])      # values, without ghosts
       x, v = np.concatenate(xs), np.concatenate(values)
       order = np.argsort(x)                   # patches may come in any order
       return x[order], v[order]


   run = Run("first_run")
   t = 50.0
   B = run.GetB(t)       # components "x", "y", "z"
   V = run.GetVi(t)      # components "x", "y", "z"
   N = run.GetNi(t)      # a single component, called "value"

   x, by = profile(B, "y")
   x, vy = profile(V, "y")
   x, n = profile(N, "value")

   fig, (ax1, ax2) = plt.subplots(nrows=2, sharex=True, figsize=(6, 5))
   ax1.plot(x, by, label="By")
   ax1.plot(x, -vy, "--", label="-Vy")
   ax1.legend()
   ax2.plot(x, n)
   ax2.set_ylabel("ion density")
   ax2.set_xlabel("x")
   fig.savefig("first_result.png")

Here the run is small and has a single patch, but the same code works for any
number of patches.

What you should see in ``first_result.png``:

* **Top:** ``-Vy`` lies on top of ``By``. In an Alfvén wave going along
  :math:`+\mathbf{B}_0`, the velocity perturbation is
  :math:`\delta\mathbf{V} = -\delta\mathbf{B}/\sqrt{\mu_0 \rho}`, which is
  :math:`-\delta\mathbf{B}` in PHARE units with density 1. You put this in the
  initial condition; the plot shows that it still holds at t = 50, after the
  wave has moved half a box.
* **Bottom:** the density stays at 1, with random fluctuations of a few
  percent (up to about 15 % at single points). An Alfvén wave does not compress the plasma; the fluctuations are
  the noise of a finite number of particles per cell.

.. admonition:: Known issue
   :class: note

   ``merged=True`` (used in :doc:`first_simulation`) currently works only with
   ``GetB`` and ``GetE``. With ``GetNi``, ``GetVi`` and the other moment getters
   it stops with ``AttributeError: 'dict' object has no attribute 'domain_box'``.
   Use the patch loop above for those quantities.


Common mistakes
---------------

* ``KeyError: 'By'`` in ``profile(B, "By")`` — without ``merged=True``, the
  components of a vector are called ``"x"``, ``"y"``, ``"z"``. With
  ``merged=True`` they are ``"Bx"``, ``"By"``, ``"Bz"``.
* ``KeyError: 'value'`` or ``'x'`` — a density has one component,
  ``"value"``; a vector has three, ``"x"``, ``"y"``, ``"z"``.
* ``FileNotFoundError: ... unable to open file: name = 'first_run/ions_pop_protons_flux.h5'``
  — that quantity was not written. Add the matching diagnostics block to your
  input script and run it again.
* ``KeyError: "Unable to synchronously open object (object '70.0000000000' doesn't exist)"``
  — nothing was written at that time (here t = 70). Use one of the times
  printed by ``run.all_times()``.


Next
----

* :doc:`/data_analysis/reading_data`: every getter, and the files they read.
* :doc:`/data_analysis/plotting_fields`: plotting in 1D and 2D.
* :doc:`/user_guide/simulation`: all the settings of the ``Simulation`` block.
