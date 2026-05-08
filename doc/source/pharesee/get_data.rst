
Loading and accessing data
==========================

PHARE simulation output is read through the ``pharesee`` package, part of ``pyphare``.
Make sure ``pyphare`` is on your Python path:

.. code-block:: bash

    export PYTHONPATH=$PYTHONPATH:/path/to/PHARE/pyphare


Creating a Run object
---------------------

The entry point for all data access is the :class:`~pyphare.pharesee.run.Run` class.
It represents a completed (or in-progress) simulation and provides getter methods for
every diagnostic quantity.

.. code-block:: python

    from pyphare.pharesee.run import Run

    run = Run("/data/scratch/phare_jobs/run001")

An optional ``default_time`` argument pins the fallback time used by helpers
such as :meth:`~pyphare.pharesee.run.Run.GetDomainSize` and
:meth:`~pyphare.pharesee.run.Run.GetDl` when no time is supplied:

.. code-block:: python

    run = Run("/data/scratch/phare_jobs/run001", default_time=10.0)

``Run`` scans the output directory on construction and raises ``RuntimeError`` if
no ``*.h5`` or ``*.vtkhdf`` files are found.


Querying available times
------------------------

Each diagnostic produces a separate file whose name matches the quantity key.
Use these methods to discover what times were written.

``run.times(key)``
    Return the array of times written for a single diagnostic.
    The *key* is the HDF5 filename stem, for example ``"EM_B"``, ``"EM_E"``,
    or ``"ions_charge_density"``.

``run.all_times()``
    Return a dict mapping every available filename stem to its array of times.

.. code-block:: python

    # times at which the magnetic field was written
    t_B = run.times("EM_B")

    # survey everything that was written
    all_t = run.all_times()
    for diag, times in all_t.items():
        print(diag, times)


Getter methods
--------------

All getters share the same positional ``time`` argument and a set of common
keyword arguments described in the :ref:`kwargs` section below.
They return a :class:`~pyphare.pharesee.hierarchy.PatchHierarchy` (or a wrapped
:class:`~pyphare.pharesee.hierarchy.ScalarField` /
:class:`~pyphare.pharesee.hierarchy.VectorField`), or — when ``merged=True`` — a
dict of interpolators keyed by component name.

Electromagnetic fields
~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Method
     - Description
   * - ``GetB(time, ...)``
     - Magnetic field **B**. Components moved to primal nodes by default
       (``all_primal=True``).
   * - ``GetE(time, ...)``
     - Electric field **E**. Same primal-node behaviour as ``GetB``.
   * - ``GetJ(time, ...)``
     - Current density **J**, computed from ``curl B``. Moved to primal nodes
       by default.
   * - ``GetDivB(time, ...)``
     - Divergence of **B** (diagnostic quantity).

Ion fluid — total
~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Method
     - Description
   * - ``GetNi(time, ...)``
     - Total ion charge density.
   * - ``GetMassDensity(time, ...)``
     - Total ion mass density.
   * - ``GetVi(time, ...)``
     - Ion bulk velocity.
   * - ``GetPi(time, ...)``
     - Ion pressure tensor (computed from momentum tensor, mass density,
       and bulk velocity).
   * - ``GetPe(time, ...)``
     - Electron pressure ``Te * Ni`` (isothermal closure; ``Te`` is read from
       the serialised simulation object).

Ion fluid — per population
~~~~~~~~~~~~~~~~~~~~~~~~~~

These methods require the population name (a string matching the name used in
the simulation configuration, e.g. ``"protons"``).

.. list-table::
   :header-rows: 1
   :widths: 45 55

   * - Method
     - Description
   * - ``GetN(time, pop_name, ...)``
     - Number density for a single ion population.
   * - ``GetFlux(time, pop_name, ...)``
     - Momentum flux for a single ion population.
   * - ``GetPressure(time, pop_name, ...)``
     - Pressure tensor for a single ion population (computed from the
       population momentum tensor, flux, and density).

MHD quantities
~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Method
     - Description
   * - ``GetMHDrho(time, ...)``
     - MHD mass density.
   * - ``GetMHDV(time, ...)``
     - MHD velocity.
   * - ``GetMHDP(time, ...)``
     - MHD pressure.
   * - ``GetMHDrhoV(time, ...)``
     - MHD momentum density.
   * - ``GetMHDEtot(time, ...)``
     - MHD total energy.

Particles
~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 45 55

   * - Method
     - Description
   * - ``GetParticles(time, pop_name)``
     - Particle data for one population (or a list/tuple of populations).
       Returns a :class:`~pyphare.pharesee.hierarchy.PatchHierarchy` whose
       patch data contain particle arrays rather than field arrays.

Metadata and diagnostics
~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Method
     - Description
   * - ``GetTags(time, ...)``
     - AMR refinement tags per patch.
   * - ``GetRanks(time, ...)``
     - MPI rank assigned to each patch (derived from the ``EM_B`` file).
   * - ``GetParticleCount(time)``
     - Per-patch particle count diagnostic.


.. _kwargs:

Common keyword arguments
------------------------

All getter methods accept the following keyword arguments.

``merged=True``
    Instead of returning a raw :class:`~pyphare.pharesee.hierarchy.PatchHierarchy`,
    project all AMR levels onto the finest uniform grid and return a dict of
    interpolators keyed by component name (e.g. ``{"Bx": interp, "By": interp,
    "Bz": interp}``).  Each interpolator is callable with coordinate arrays and
    returns the field values on those coordinates.

``interp="nearest"`` or ``interp="linear"``
    Interpolation method used when projecting levels during merging (only
    relevant when ``merged=True``).  Defaults to ``"nearest"``.

``all_primal=True``
    Move staggered (Yee-grid) data to primal (cell-corner) nodes before
    returning.  Applies to ``GetB``, ``GetE``, ``GetJ``, and MHD getters.
    Automatically disabled when ``merged=True`` (merging works on the native
    staggered layout).

``selection_box``
    Restrict loading to a spatial sub-region.  Pass a four-element tuple
    ``(x_lo, y_lo, x_hi, y_hi)`` (2-D) or the appropriate lower/upper
    coordinates for your dimensionality.  Only patches that overlap the box
    are loaded from disk.

.. code-block:: python

    # Load only a sub-region of the domain
    B = run.GetB(10.0, selection_box=(20.0, 0.0, 40.0, 10.0))


The PatchHierarchy
------------------

Getter methods return a :class:`~pyphare.pharesee.hierarchy.PatchHierarchy`
that organises data as a tree of levels, each containing a list of patches.

Navigating levels and patches
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    hier = run.GetB(10.0)

    # dict: level index -> PatchLevel
    levels = hier.levels(10.0)

    for ilvl, level in levels.items():
        print(f"Level {ilvl}: {len(level.patches)} patches")
        for patch in level.patches:
            print(f"  patch box: {patch.box}")

Accessing field data on a patch
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Each :class:`~pyphare.pharesee.hierarchy.patch.Patch` exposes its data through
the ``patch_datas`` dict, keyed by component name:

.. code-block:: python

    for ilvl, level in hier.levels(10.0).items():
        for patch in level.patches:
            pd = patch.patch_datas["By"]
            x = pd.x          # 1-D coordinate array (primal or dual)
            data = pd.dataset  # numpy array of field values

Listing available quantities
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    hier.quantities()   # e.g. ["Bx", "By", "Bz"]

Complete 1-D extraction example
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    import numpy as np

    hier = run.GetB(10.0, all_primal=True)

    all_x, all_By = [], []
    for patch in hier.level(0).patches:
        pd = patch.patch_datas["By"]
        all_x.append(pd.x)
        all_By.append(pd.dataset[:])

    x  = np.concatenate(all_x)
    By = np.concatenate(all_By)
    idx = np.argsort(x)
    x, By = x[idx], By[idx]


Working with merged data
------------------------

When ``merged=True`` the getter returns a dict of interpolators instead of a
hierarchy.  The interpolators can be evaluated at any set of coordinates within
the domain.

.. code-block:: python

    import numpy as np
    import matplotlib.pyplot as plt

    # Project all levels of B onto the finest grid
    merged_B = run.GetB(10.0, merged=True, interp="linear")

    # merged_B is {"Bx": (interp_fn, coords), "By": ..., "Bz": ...}
    By_interp, (x,) = merged_B["By"]   # 1-D: coords is a 1-tuple

    plt.plot(x, By_interp(x))
    plt.xlabel("x")
    plt.ylabel("By")
    plt.show()

For 2-D runs, ``coords`` is a 2-tuple ``(x, y)`` and the interpolator accepts
2-D coordinate arrays:

.. code-block:: python

    By_interp, (x, y) = merged_B["By"]
    X, Y = np.meshgrid(x, y, indexing="ij")
    By_2d = By_interp(X, Y)


Direct file loading
-------------------

For one-off loading without constructing a ``Run`` object, use the
:func:`~pyphare.pharesee.hierarchy.hierarchy_from` function directly.  It
dispatches automatically based on file extension (``*.h5`` or ``*.vtkhdf``):

.. code-block:: python

    from pyphare.pharesee.hierarchy import hierarchy_from

    hier = hierarchy_from(h5_filename="EM_B.h5", times=[10.0])

Pass ``hier=existing_hier`` to accumulate multiple quantities into a single
hierarchy object.


API reference
-------------

.. autoclass:: pyphare.pharesee.run.Run
   :members:
   :undoc-members:
   :show-inheritance:
