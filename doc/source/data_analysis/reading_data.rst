================
Reading the data
================

**Goal:** read the files written by a PHARE run back into Python, with the
``Run`` object of the ``pharesee`` package.

The examples use the 1D run of :doc:`/getting_started/first_simulation`, with
more diagnostics added (see `Minimal example`_).

.. contents:: On this page
   :local:
   :depth: 1


Setting up
----------

The data produced by PHARE is read with the ``pharesee`` Python package, part
of ``pyphare``. Make sure your ``PYTHONPATH`` includes the folder where
``pyphare`` is located.

.. code-block:: bash

    export PYTHONPATH=$PYTHONPATH:/path/to/PHARE/pyphare

Importing ``pharesee`` should then work:

.. code-block:: python

    import pyphare.pharesee as pharesee


Minimal example
---------------

Take ``first_simulation.py`` from :doc:`/getting_started/first_simulation` and
replace its single ``ElectromagDiagnostics`` line by these lines (each
``...Diagnostics`` line asks PHARE to write one quantity; see
:doc:`/user_guide/diagnostics`):

.. code-block:: python

   times = [0.0, 25.0, 50.0]
   ph.ElectromagDiagnostics(quantity="B", write_timestamps=times)
   ph.ElectromagDiagnostics(quantity="E", write_timestamps=times)
   ph.FluidDiagnostics(quantity="charge_density", write_timestamps=times)
   ph.FluidDiagnostics(quantity="bulkVelocity", write_timestamps=times)
   ph.FluidDiagnostics(quantity="density", population_name="protons",
                       write_timestamps=times)
   ph.FluidDiagnostics(quantity="flux", population_name="protons",
                       write_timestamps=times)

Run it, then read the results. The first step is to create a ``Run`` object. It
represents your simulation and lets you reach its data. Give it the folder
that the simulation wrote (the ``dir`` of ``diag_options``):

.. code-block:: python

    from pyphare.pharesee.run import Run

    run = Run("first_run")
    print(run.times("B"))        # [ 0. 25. 50.]
    print(run.GetDomainSize())   # [100.]
    print(run.GetDl())           # [0.5]  cell size on level 0

    B = run.GetB(25.0)           # a patch hierarchy (see below)
    Ni = run.GetNi(25.0)         # ion density, also a patch hierarchy

Every data getter has the same form, ``run.GetXXX(time)``, where ``XXX`` is the
name of a physical quantity and ``time`` is the time you want. The time must be
one at which the quantity was written: any other time gives a ``KeyError``
(``Unable to synchronously open object (object '7.0000000000' doesn't
exist)``). Use ``run.times("B")`` to list the available times of a quantity.
``run.all_times()`` gives them for every file, in a dictionary.


What it does
------------

A getter opens the one file that holds the quantity (written by a
:doc:`diagnostics </user_guide/diagnostics>` block), reads the requested time,
and returns a Python object: a :doc:`patch hierarchy <patch_hierarchy>` by
default, or an interpolator and its grid with ``merged=True``.

Which getter reads which file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The file name is fixed by the diagnostic you asked for. The getter fails with
``FileNotFoundError`` if the diagnostic was not requested in the input script.

.. list-table::
   :header-rows: 1
   :widths: 24 38 38

   * - Getter
     - Diagnostic to ask for
     - File read
   * - ``GetB(t)``
     - ``ElectromagDiagnostics("B")``
     - ``EM_B.h5``
   * - ``GetE(t)``
     - ``ElectromagDiagnostics("E")``
     - ``EM_E.h5``
   * - ``GetNi(t)``
     - ``FluidDiagnostics("charge_density")``
     - ``ions_charge_density.h5``
   * - ``GetMassDensity(t)``
     - ``FluidDiagnostics("mass_density")``
     - ``ions_mass_density.h5``
   * - ``GetVi(t)``
     - ``FluidDiagnostics("bulkVelocity")``
     - ``ions_bulkVelocity.h5``
   * - ``GetChargeFlux(t)``
     - ``FluidDiagnostics("charge_flux")``
     - ``ions_charge_flux.h5``
   * - ``GetN(t, pop)``
     - ``FluidDiagnostics("density", population_name=pop)``
     - ``ions_pop_<pop>_density.h5``
   * - ``GetFlux(t, pop)``
     - ``FluidDiagnostics("flux", population_name=pop)``
     - ``ions_pop_<pop>_flux.h5``
   * - ``GetPressure(t, pop)``
     - per-population ``"momentum_tensor"``, ``"flux"`` and ``"density"``
     - ``ions_pop_<pop>_momentum_tensor.h5`` + the two above
   * - ``GetPi(t)``
     - ``"momentum_tensor"``, ``"mass_density"`` and ``"bulkVelocity"`` (ions)
     - ``ions_momentum_tensor.h5`` + the two above
   * - ``GetPe(t)``
     - ``FluidDiagnostics("charge_density")``; result is ``Te`` times the density
     - ``ions_charge_density.h5``
   * - ``GetJ(t)``
     - ``ElectromagDiagnostics("B")``; current computed from B
     - ``EM_B.h5``
   * - ``GetDivB(t)``
     - ``ElectromagDiagnostics("B")``; computed from B
     - ``EM_B.h5``
   * - ``GetParticles(t, pop)``
     - ``ParticleDiagnostics("domain", population_name=pop)``, see :doc:`distributions`
     - ``ions_pop_<pop>_domain.h5``
   * - ``GetParticleCount(t)``
     - ``InfoDiagnostics("particle_count")``
     - ``particle_count.h5``
   * - ``GetMHDrho``, ``GetMHDV``, ``GetMHDP``, ``GetMHDrhoV``, ``GetMHDEtot``
     - ``MHDDiagnostics`` with quantity ``rho``, ``V``, ``P``, ``rhoV``, ``Etot``
     - ``mhd_rho.h5``, ``mhd_V.h5``, ``mhd_P.h5``, ``mhd_rhoV.h5``, ``mhd_Etot.h5``

The other methods do not read a field: ``GetDomainSize()`` (box size),
``GetDl(level=0)`` (cell size on a level), ``GetMass(pop)`` and
``GetCharge(pop)`` (read from the particle population attributes in the files),
``times(name)`` and ``all_times()``.

The names accepted by ``times(name)`` are short names for some files
(``"B"``, ``"E"``, ``"Ni"``, ``"Vi"``, ``"rho"`` for the mass density, ``"CF"``
for the charge flux, ``"nppc"``, ``"mhd_rho"``...); for the other files it is
the file name without ``.h5``, for example
``run.times("ions_pop_protons_flux")``.

Two kinds of results
^^^^^^^^^^^^^^^^^^^^

**By default** a getter returns a *patch hierarchy*: PHARE stores its data in
pieces of grid called patches, grouped by refinement level, and the hierarchy
keeps that structure. It is what you need to plot with
``hierarchy.plot(...)`` (see :doc:`plotting_fields`) and to work on several
levels. Its structure is explained in :doc:`/data_analysis/patch_hierarchy`.

.. code-block:: python

    B = run.GetB(25.0)       # B is a patch hierarchy
    print(B.quantities())    # ['x', 'y', 'z']

**With** ``merged=True`` the patches are put back together on the finest grid
available, and you get, for each component, a function you can evaluate and the
grid points where it is defined:

.. code-block:: python

    Bm = run.GetB(25.0, merged=True)   # a dictionary: "Bx", "By", "Bz"
    by, (x,) = Bm["By"]                # by is a function, x the grid points
    values = by(x)

In 1D, ``x`` is one array. In 2D, ``(x, y)`` are two 1D arrays, and the function
takes 2D arrays: use ``X, Y = numpy.meshgrid(x, y, indexing="ij")`` then
``by(X, Y)``.

.. note:: **Known issue.** ``merged=True`` currently works only for ``GetB`` and
   ``GetE``. For ion quantities (``GetNi``, ``GetVi``, ``GetN``, ``GetFlux``...)
   it stops with ``AttributeError: 'dict' object has no attribute 'domain_box'``
   (or a ``KeyError``): use the patch hierarchy instead, as in
   :doc:`plotting_fields`. ``GetPe`` ignores ``merged=True`` (it returns a
   hierarchy). The last point of the merged grid lies half a cell past the box
   (`#1553 <https://github.com/PHAREHUB/PHARE/issues/1553>`_, merged coordinates
   overshoot), so its value is extrapolated.


Options
-------

These options are shared by most getters.

.. list-table::
   :header-rows: 1
   :widths: 22 56 22

   * - Option
     - Meaning
     - Default
   * - ``time``
     - time to read (first argument); must be a written time
     - none
   * - ``merged``
     - ``True``: functions on the finest grid (see the note above); ``False``:
       patch hierarchy
     - ``False``
   * - ``interp``
     - interpolation used by ``merged=True``. 1D: any ``kind`` of
       ``scipy.interpolate.interp1d`` (for example ``"nearest"``, ``"linear"``);
       2D: ``"nearest"`` or ``"bilinear"``
     - ``"nearest"``
   * - ``all_primal``
     - for ``GetB``, ``GetE``, ``GetJ`` and the MHD getters: move all
       components to the same grid points (primal nodes, the cell corners).
       Components of B and E live at different points inside a cell; this puts
       them together. Forced to ``False`` when ``merged=True``
     - ``True``

``Run(path)`` takes the folder of the run (its ``default_time`` argument is
currently unused).


Common mistakes
---------------

* ``FileNotFoundError: ... name = 'first_run/ions_pop_alpha_domain.h5' ...`` — the
  population name is wrong (here there is no ``alpha``), or the diagnostic was
  not requested. Use the population names of your input script.
* ``FileNotFoundError: ... 'first_run/ions_pop_protons_momentum_tensor.h5' ...``
  — ``GetPressure`` needs a per-population ``momentum_tensor`` diagnostic.
* ``KeyError: "Unable to synchronously open object (object '7.0000000000'
  doesn't exist)"`` — time 7.0 was not written. List the times with
  ``run.times("B")``.
* ``ValueError: divB is 0 by construction in 1D`` — ``GetDivB`` only makes
  sense in 2D and 3D.
* In 1D, ``GetJ(t)`` stops with ``KeyError: 'Jx'``, and ``GetPressure`` and
  ``GetPi`` with an ``AssertionError``. They work in 2D. ``GetJ(t, all_primal=False)``
  works in 1D.
* ``Run`` only looks at ``*.h5`` files in the folder. Use the ``phareh5``
  format for analysis with ``pharesee``.

.. note:: Diagnostics are written in single precision, so fields carry about
   7 significant digits, unless PHARE was built with ``-DPHARE_DIAG_DOUBLES=1``
   (see :doc:`/getting_started/build`).

Next: :doc:`plotting_fields` and :doc:`distributions`.


Reference
---------

.. autoclass:: pyphare.pharesee.run.Run
    :members:
    :undoc-members:
    :show-inheritance:
