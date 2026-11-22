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
of ``pyphare``. Use the same ``PYTHONPATH`` line as for running (see
:doc:`/getting_started/build`):

.. code-block:: bash

    export PYTHONPATH=/path/to/PHARE/pyphare:/path/to/build:$PYTHONPATH

Reading data needs only the first folder, so you can also analyse a run on a
machine where PHARE is not built.

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
one at which the quantity was written: any other time raises a
``TimeNotFoundError`` (a kind of ``KeyError``) whose message lists the nearest
written times. Use ``run.times("B")`` to list the available times of a
quantity, sorted. ``run.all_times()`` gives them for every file, in a
dictionary.


What it does
------------

A getter opens the one file that holds the quantity (written by a
:doc:`diagnostics </user_guide/diagnostics>` block), reads the requested time,
and returns a Python object: a :doc:`patch hierarchy <patch_hierarchy>` by
default, or an interpolator and its grid with ``merged=True``.

The default getters go through the AMR reader (:doc:`amr_reader`): the
hierarchy reads its patches only when you use them, and a getter also accepts
``box=`` to read part of the domain. ``merged=True``, ``all_primal=False``,
particles and pressure tensors are still read by the previous reader; the
:doc:`amr_reader` page lists what differs between the two.

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
   * - ``GetB0(t)``
     - ``ElectromagDiagnostics("B0")``, with an
       :doc:`external field </user_guide/external_field>`
     - ``EM_B0.h5``
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
(``"B"``, ``"E"``, ``"B0"``, ``"Ni"``, ``"Vi"``, ``"rho"`` for the mass density, ``"CF"``
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

.. note:: **Known issue.** ``merged=True`` currently works only for ``GetB``,
   ``GetE`` and ``GetB0``. For ion quantities (``GetNi``, ``GetVi``, ``GetN``, ``GetFlux``...)
   it stops with ``AttributeError: 'dict' object has no attribute 'domain_box'``
   (or a ``KeyError``): use the patch hierarchy instead, as in
   :doc:`plotting_fields`. ``GetPe`` ignores ``merged=True`` (it returns a
   hierarchy). The last point of the merged grid lies half a cell past the box
   (`#1553 <https://github.com/PHAREHUB/PHARE/issues/1553>`_, merged coordinates
   overshoot), so its value is extrapolated.


Reconnection rate (2D)
----------------------

Three methods of ``Run`` measure magnetic reconnection in a 2D run. They read
``EM_B.h5`` (``ElectromagDiagnostics("B")``) and work on the finest level,
interpolated over the whole box.

.. code-block:: python

    import numpy as np

    times = np.sort(run.times("B"))          # sorted NumPy array, not a list
    tc, rate, flux, xpoint = run.GetReconnectionRate(times)

* ``run.GetMagneticFlux(t)`` returns ``Az, (x, y)``: the flux function
  :math:`A_z`, with :math:`B_x = \partial A_z/\partial y` and
  :math:`B_y = -\partial A_z/\partial x`, on the nodes of the finest grid
  (from 0 to the box size included), with :math:`A_z = 0` at the origin.
* ``run.FindPrimaryXPoint(Az, x, y)`` returns ``x_X, y_X, index``: among the
  5 % of nodes where :math:`|\nabla A_z|` is smallest, the saddle point of
  :math:`A_z` with the most negative Hessian determinant.
* ``run.GetReconnectionRate(times)`` does both at each time and returns the
  mid-times ``tc`` (one fewer than ``times``), the rate
  :math:`\mathrm{d}A_z/\mathrm{d}t` at the X point between two outputs, the
  flux ``flux`` at the X point at each time, and the X-point positions
  ``xpoint`` (an array of shape ``(len(times), 2)``).

The rate is in code units, :math:`B_0 v_A`, that is the reconnection electric
field normalized to the field and Alfvén speed of the normalization (see
:doc:`/theory/normalization`); its sign follows the orientation of
:math:`A_z`.

.. warning::

   The X point is searched over the whole box, at each time separately. With
   several X points, for example the two current sheets of a periodic Harris
   setup, it can jump from one to another between outputs, and the rate then
   mixes them. Check ``xpoint`` before trusting ``rate``. These methods are not
   covered by PHARE's tests.

``times`` can be a list or a NumPy array. With the previous reader
(``Run(path, legacy=True)``) it must be a NumPy array: a Python list stops with
``TypeError: unsupported operand type(s) for /: 'list' and 'int'``.


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
       2D: ``"nearest"`` or ``"bilinear"``. Without ``merged=True`` it raises
       ``TypeError``
     - ``"nearest"``
   * - ``box``
     - read only the patches meeting this box, in physical units:
       ``box=dict(x=(10.0, 30.0))`` (add ``y=``, ``z=`` in 2D and 3D). Not with
       ``merged=True`` or ``all_primal=False``
     - the whole domain
   * - ``all_primal``
     - for ``GetB``, ``GetE``, ``GetB0``, ``GetJ`` and the MHD getters: move all
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
* ``TimeNotFoundError: first_run/EM_B.h5: no data at time 7.0; nearest times:
  0.0, 25.0, 50.0 ...`` — time 7.0 was not written. List the times with
  ``run.times("B")``.
* In 1D, ``GetDivB`` returns zeros, and NaN next to a coarse-fine boundary:
  divB is 0 by construction (the previous
  reader, ``legacy=True``, stops with ``ValueError``).
* In 1D, ``GetPressure`` and ``GetPi`` stop with an ``AssertionError``. They
  work in 2D.
* ``TypeError`` naming an unexpected keyword — a typo in a getter keyword; the
  message lists the accepted ones.
* ``Run`` reads the ``*.h5`` (``phareh5``) and ``*.vtkhdf`` (``pharevtkhdf``)
  files of the folder. Particles, pressure tensors and ``merged=True`` need
  ``phareh5`` files.

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
