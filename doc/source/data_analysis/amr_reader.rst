==============
The AMR reader
==============

**Goal:** read a large AMR run in seconds, from ``.h5`` or ``.vtkhdf`` files, with
the calls of :doc:`reading_data`; know what differs from the previous reader and
how to go back to it.

.. contents:: On this page
   :local:
   :depth: 1


Minimal example
---------------

The getters of :doc:`reading_data` keep their names and arguments. These calls
are new (on a 2D run):

.. code-block:: python

    from pyphare.pharesee.run import Run

    run = Run(path)             # a folder of .h5 or .vtkhdf files
    print(run)                  # says which reader is in use
    t = 0.1

    # read only the patches meeting a box (physical units); the box becomes the field's extent
    B = run.GetB(t, box=dict(x=(5.0, 15.0), y=(3.0, 8.0)))

    B.z.plot()                  # one image, same keywords as before (+ robust=, max_level=, ...)
    grid = B.z.finest()         # whole box on the finest node grid
    grid.data, grid.level_map   # values [x, y] and the AMR level each node shows

    J = run.GetJ(t)             # computed on demand, also takes box=
    x, bz = B.z(y=5.0)          # line cut at y = 5
    B.z.mean(), B.z.integral()  # reductions that count each point once, finest data
    run.plan("B", t)            # what a read would touch (patches, bytes), without reading

To go back to the previous reader:

.. code-block:: python

    run = Run(path, legacy=True)     # one run

.. code-block:: bash

    export PHARESEE_LEGACY=1         # every Run in the process

This works only for ``.h5`` files; ``.vtkhdf`` files raise a ``ValueError``.


What it does
------------

``Run`` reads ``.h5`` and ``.vtkhdf`` diagnostics with a reader that keeps each
AMR level in flat arrays and reads only the patches a request needs
(``pyphare.pharesee.amr``). It is the default for ``.h5`` files and the only
reader for ``.vtkhdf`` files. A getter returns a field that reads its data when
you first use it: creating ``B = run.GetB(t)`` costs almost nothing.

What the new reader cannot serve yet (``merged=True``, ``all_primal=False``,
particles, pressure tensors, tags, ranks) is handed to the previous reader
automatically, so those calls return the previous reader's objects.
:doc:`/dev/pharesee_reader` lists which reader serves each call. The previous
reader stays only until particles, pressure tensors, tags and ranks are ported.

Differences from the previous reader
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Same calls, but these behaviours changed on purpose. Tell us if one blocks you.

.. list-table::
   :header-rows: 1
   :widths: 26 26 48

   * - What
     - Previous reader
     - New reader (why)
   * - Times (``B.times()``, ``B.time_hier`` keys)
     - text ``"0.0500000000"``
     - numbers; ``time_hier["0.0500000000"]`` still works (numbers sort and compare)
   * - ``run.times("B")``
     - in the order of the file
     - sorted
   * - Time not written (``GetB(7.0)``)
     - ``TimeNotFoundError``
     - the same class, its message lists the nearest written times
   * - Precision
     - float64
     - the file's (float32 unless PHARE was built with ``PHARE_DIAG_DOUBLES``): half the
       memory, the previous values rounded to it (see :ref:`dev-pharesee_reader-numerics`)
   * - Ghost nodes of per-patch ``dataset``
     - file values
     - ``.h5``: the same; ``.vtkhdf``: NaN (the file has no ghosts)
   * - Per-patch ``dataset``
     - writeable
     - read-only: the data is shared between fields; ``copy.deepcopy(patch_data)`` or
       ``field.copy()`` give writeable data
   * - Unknown keyword (``GetB(t, foo=1)``)
     - ignored
     - ``TypeError`` (a typo no longer passes silently)
   * - ``interp=`` without ``merged=True``
     - ignored
     - ``TypeError`` (``interp`` only applies to ``merged=True``)
   * - ``GetB([])``
     - ``IndexError``
     - ``ValueError`` saying a time is needed
   * - ``field == 1.0``, ``bool(field)``, ``iter(field)``
     - a single bool / empty iteration
     - ``TypeError``: one truth value of AMR data (overlapping levels) is ambiguous;
       compare with ``np.isclose``/``np.equal``, then ``np.all``/``np.any``; or
       ``field.max()``, ``field.times()``
   * - ``B.levels()`` or ``B.plot()`` on a field of several times
     - the first time, silently
     - ``ValueError`` asking for ``time=``
   * - 1D ``plot``
     - one line and one legend entry per patch; returns nothing
     - one line per level (``L0``, ``L1``); returns ``(fig, ax)`` as in 2D
   * - ``plot(qty="By")`` on ``GetB(t)``
     - ``KeyError`` (the components are ``"x"``, ``"y"``, ``"z"``)
     - accepted, as ``qty="y"``
   * - 3D ``plot``
     - draws nothing
     - ``ValueError``: take a cut first (``finest(box=...)`` or a line cut)
   * - ``GetDivB`` in 1D
     - ``ValueError``
     - zeros (NaN next to a coarse-fine boundary)
   * - ``GetJ`` in 1D
     - ``KeyError: 'Jx'``
     - works
   * - Line cut ``B("x", y=5.0)``
     - tuple
     - ``B.x(y=5.0)``
   * - ``hierarchy_compare(new, previous)``
     -
     - data mismatch at the last float32 digit (the precision row above); compare with
       ``np.allclose`` at float32 tolerance
   * - ``VectorField.FROM``, ``dist_plot``, ``resolve_selection_boxes``
     - available
     - not available (``isinstance(B, PatchHierarchy)`` still holds while both readers
       coexist)
   * - ``selection_box=``
     - integer box
     - works, with a ``DeprecationWarning``: prefer ``box=dict(x=(x0, x1), ...)``
       (physical units, float bounds)
   * - ``GetMagneticFlux``
     - merged finest data
     - finest data painted with ``finest()``: within ~1-2 % of ``max|Az|``


Options
-------

Getters
^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 30 50 20

   * - Option
     - Meaning
     - Default
   * - ``box``
     - every field getter: read only the patches meeting this box, in physical units,
       ``dict(x=(x0, x1), y=(y0, y1))``; the box becomes the field's extent
     - the whole domain
   * - ``stencil`` (``GetJ``)
     - ``"yee"``: the code's J, from the Yee B (``.h5`` only); ``"node"``: the curl of
       node B
     - ``"yee"`` on ``.h5``, ``"node"`` on ``.vtkhdf``
   * - ``stencil`` (``GetDivB``)
     - ``"yee"``: the previous reader's per-cell divB (2D ``.h5``); ``"node"``: at
       the nodes
     - ``"yee"`` on 2D ``.h5``, else ``"node"``
   * - ``mask_faces`` (``GetDivB``, ``stencil="node"``)
     - NaN at nodes next to a coarse-fine face, where a one-sided difference or ghost
       values were needed
     - ``True``

Plot and finest grid
^^^^^^^^^^^^^^^^^^^^

``plot`` takes the keywords of :doc:`plotting_fields`, plus:

.. list-table::
   :header-rows: 1
   :widths: 30 50 20

   * - Option
     - Meaning
     - Default
   * - ``robust``
     - colour limits from the 1st and 99.7th percentiles
     - ``False``
   * - ``symmetric``
     - colour limits symmetric around 0
     - ``False``
   * - ``norm``
     - a matplotlib ``Normalize``, instead of ``vmin``/``vmax``
     - none
   * - ``max_level``
     - finest level drawn
     - matched to the axes (below)
   * - ``fill``
     - how nodes no level holds are filled, as in ``finest(fill=...)`` below
     - ``"cubic"``
   * - ``interpolation``, ``interpolation_stage``
     - passed to matplotlib's ``imshow``
     - ``"antialiased"``

A 2D ``plot()`` builds its picture at a resolution matched to the axes (at most
2 nodes per pixel), so a whole-domain panel does not alias thin structures and
does not read the finest level everywhere. Zoom with ``xlim=``/``ylim=`` to get
the finest level back, or force a level with ``max_level=``. The default colour
limits span every level in view, finite values only.

``finest(fill=...)`` fills the fine nodes a finer level does not cover from the
coarser level: ``"cubic"`` (default) is a Catmull-Rom interpolation clipped to
its two neighbours; ``"nearest"`` copies a coarse node. These values are for
display only: reductions (``mean``, ``integral``, ``min``, ``max``) use the data.

Environment variables
^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 28 28 44

   * - Variable
     - Default
     - Effect
   * - ``PHARESEE_LEGACY``
     - unset (new reader)
     - ``1/true/yes``: previous reader for ``.h5``; ``0/false/no``: new reader
   * - ``PHARESEE_IO_THREADS``
     - 8
     - Reader threads (at least 1)
   * - ``PHARESEE_CACHE``
     - on
     - ``0/false/no/off`` disables the on-disk cache of ``.h5`` indexes
   * - ``PHARESEE_CACHE_DIR``
     - ``$XDG_CACHE_HOME/pharesee``, else ``~/.cache/pharesee``
     - Directory of that cache
   * - ``PHARESEE_CACHE_MB``
     - 1024
     - Size cap (MiB) of that cache
   * - ``PHARESEE_GEOM_CACHE_MB``
     - 1024
     - Size cap (MiB) of the in-memory cache of level layouts and face tables
   * - ``PHARESEE_RAW_IO``
     - on
     - ``0`` reads through h5py instead of direct reads: same results, slower
   * - ``PHARESEE_SAMPLE_DATA``
     - ``<repo>/build/sample_data``
     - Tests only: where the sample runs are


Common mistakes
---------------

* ``TypeError: GetB(interp='linear'): interp only applies to merged=True ...`` —
  ``interp`` does nothing without ``merged=True``: drop it.
* ``TypeError`` naming an unexpected keyword — a typo in a getter or ``plot``
  keyword; the message lists the accepted ones.
* ``ValueError: VectorField EM_B holds 2 times [0.0, 0.05]; pass time=...`` — a
  field read at several times (``run.GetB([t0, t1])``) needs ``time=`` in
  ``plot``, ``levels`` and the reductions, or ``B.at(t)``.
* ``ValueError: assignment destination is read-only`` (writing into a
  per-patch ``dataset``) or ``AttributeError: the dataset of a file-backed patch
  view is read-only`` (replacing it) — copy first: ``copy.deepcopy(patch_data)``
  or ``field.copy()``.
* ``TypeError`` from ``if B.x == 1.0:`` or ``if B:`` — compare element by element
  with ``np.isclose(B.x, 1.0)`` or ``np.equal(B.x, 1.0)``, then reduce with
  ``np.all``/``np.any``; or use ``B.x.max()``.
* ``ValueError`` on a ``.vtkhdf`` run with ``legacy=True``, ``merged=True``,
  ``all_primal=False``, ``GetPe``, ``GetPressure``, ``GetParticles`` — the
  previous reader cannot read vtkhdf, and these quantities are not in the file;
  write ``phareh5`` diagnostics for them.
* ``mean()`` or ``integral()`` of ``GetDivB`` is NaN — the masked coarse-fine
  faces are NaN (``min``/``max`` ignore them); pass ``mask_faces=False``.

Next: how the reader decides and computes, in :doc:`/dev/pharesee_reader`.


Reference
---------

.. autoclass:: pyphare.pharesee.amr.field.ScalarField
    :members: plot, finest, mean, integral, min, max, at, times, copy

.. autoclass:: pyphare.pharesee.amr.field.VectorField
    :members: plot, at, times, copy
