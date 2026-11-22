.. _dev-pharesee_reader:

=====================
The AMR reader inside
=====================

.. container:: breadcrumb

   You are here: :doc:`/dev/index` › The AMR reader. User view:
   :doc:`/data_analysis/amr_reader`.

**In one sentence:** ``Run`` (``pyphare/pyphare/pharesee/run/run.py``) sends each
getter call either to the AMR reader (``pyphare/pyphare/pharesee/amr``) or to the
previous per-patch reader (``pyphare/pyphare/pharesee/hierarchy``), and the AMR
reader's defaults return the previous reader's values, rounded to the file's
precision.

Which reader serves what
------------------------

``Run`` decides per call. "previous" means an object of the previous reader (a
``PatchHierarchy`` of the ``hierarchy`` classes); "new" means a field of the AMR
reader.

.. list-table::
   :header-rows: 1
   :widths: 40 30 30

   * - Call
     - ``.h5``
     - ``.vtkhdf``
   * - ``GetB``, ``GetE``, ``GetB0``, ``GetNi``, ``GetN``, ``GetVi``, ``GetFlux``,
       ``GetChargeFlux``, ``GetMassDensity``, ``GetMHD*``
     - new
     - new
   * - the same with ``merged=True``
     - previous
     - ``ValueError``
   * - the same with ``all_primal=False``
     - previous (Yee-centred)
     - ``ValueError`` (the file holds nodes only)
   * - ``GetJ``
     - new; ``stencil="yee"`` (default, the code's J) or ``"node"``;
       ``all_primal=False``: previous
     - new; ``"node"`` only
   * - ``GetDivB``
     - 2D: ``stencil="yee"`` (default) gives the previous reader's per-cell divB;
       ``stencil="node"`` and 1D/3D: new
     - new, ``"node"``
   * - ``GetPe``
     - new; ``all_primal=False``: previous
     - ``ValueError`` (no simulation parameters in the file)
   * - ``GetMagneticFlux``, ``GetReconnectionRate``
     - new
     - new
   * - ``GetPressure``, ``GetPi``, ``GetTags``, ``GetRanks``, ``GetParticles``,
       ``GetParticleCount``
     - previous
     - ``ValueError`` (not written in vtkhdf)
   * - anything, with ``Run(path, legacy=True)``
     - previous
     - ``ValueError`` (the previous reader cannot read vtkhdf)
   * - anything, with ``PHARESEE_LEGACY=1``
     - previous
     - new (the variable applies to ``.h5`` only)

The node divergence (``stencil="node"``) has ``mask_faces=True`` by default: nodes
next to a coarse-fine face, where a one-sided difference or ghost values were
needed, are NaN (see `Numerics`_).

.. _dev-pharesee_reader-numerics:

Numerics
--------

* ``.h5`` files: the defaults return the previous reader's values. ``GetJ`` is the Yee
  ``J`` moved to nodes (the previous reader's values rounded to the file's precision),
  ``GetDivB`` (2D) is the previous reader's per-cell divergence, ``GetMagneticFlux`` is
  computed on the finest level's grid (within a few % of ``max|Az|``: coarse-fine gaps are
  filled differently).
* Node values averaged from Yee or cell-centred data (B, E, MHD quantities): on float32
  files, the exact mean of the cells rounded once to float32. On float64 files
  (``PHARE_DIAG_DOUBLES``), both readers round in float64 but add the cells in a different
  order, so they can differ in the last bits; where the mean cancels to ~0 that difference
  is large relative to the value, never relative to the cells.
  ``GetJ(t, stencil="node")`` gives the node stencil described below.
* ``.vtkhdf`` files only hold B at nodes. ``GetJ`` is therefore the curl of node B,
  which equals the Yee ``J`` of the code smoothed by a (1, 2, 1)/4 filter along the
  derivative direction: each wavenumber is multiplied by :math:`\cos^2(k\,\Delta x / 2)`.
  A current sheet one cell thick loses about 30 % of its peak.
* ``GetDivB`` on ``.vtkhdf`` is masked (NaN) at coarse-fine faces by default, where a node
  file has no data on one side; ``mask_faces=False`` gives values everywhere. Because of
  the NaN, ``mean()`` and ``integral()`` of a masked field are NaN; ``min()`` and
  ``max()`` ignore it.

Related
-------

* :doc:`/data_analysis/reading_data`, :doc:`/data_analysis/amr_reader` — the user view.
* :doc:`/dev/diagnostics_restarts` — how the files are written.
