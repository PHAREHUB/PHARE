.. _dev-how_to_add:

==========
How to add
==========

.. container:: breadcrumb

   You are here: :doc:`/dev/index` › :doc:`/dev/working_on_phare` › How to add.

**In one sentence:** the files to touch, in order, to add a boundary condition, a diagnostic or a
hybrid field, with the test to copy. A skipped step often fails only at run time, with the message
quoted here.

Add a boundary condition
------------------------

.. flow::
   :horizontal:

   Simulation(domain_boundaries=...) | pharein/boundary.py
   -
   dict | simulation/grid/domain_boundaries/<loc>/type
   -
   BoundaryFactory::create | one recipe per type
   -
   field refine patch strategy | fills physical ghost cells | /dev/ghost_cells

* A **boundary type** is what users write: ``"open"``, ``"reflective"``,
  ``"super-magnetofast-inflow"``. ``"none"`` is internal, for undeclared locations
  (:src:`pyphare/pyphare/pharein/boundary.py#L194-L215`). A type is a *recipe*: for each quantity, which
  **field condition** (``Dirichlet``, ``Neumann``, ``Symmetric``, ...) fills the
  :term:`ghost cells <ghost cell>`.
* MHD only at this commit: Python refuses physical boundaries for hybrid
  (:src:`pyphare/pyphare/pharein/boundary.py#L187-L192`); the hybrid model's ``BoundaryManager`` has no field
  quantities (:src:`src/amr/physical_models/hybrid_model.hpp#L94-L103`); no particle boundary yet.

New boundary type:

#. Add it at the end of ``enum class BoundaryType`` (:src:`src/core/boundary/boundary_defs.hpp#L16`).
#. Expose it to Python: ``.value("my_name", ...)`` in :src:`src/python3/cpp_etc.cpp#L165-L169`.
#. A ``Boundary`` dataclass with ``type = "my-name"`` in ``pharein/boundary.py``, declared *above*
   ``_type_to_class`` (built from the classes defined before it,
   :src:`pyphare/pyphare/pharein/boundary.py#L150-L156`).
#. A ``case`` and a ``register_my_name_conditions_`` in ``BoundaryFactory::create``
   (:src:`src/core/boundary/boundary_factory.hpp#L76-L87`), registering a condition for *every* quantity:
   a missing one stops the run with ``Field boundary condition not found.``
#. Test: a ``checkDispatch`` case in
   ``tests/core/numerics/boundary_condition/test_field_boundary_conditions_boundary_factory.cpp``.

New field condition (only if none fits): an entry in ``FieldBoundaryConditionType``, a
``field_boundary_condition_<kind>.hpp`` deriving ``IFieldBoundaryCondition``, a
``FieldBoundaryConditionSelector`` specialisation
(:src:`src/core/numerics/boundary_condition/field_boundary_condition_resolver.hpp#L18-L75`), a test copied
from ``test_field_boundary_conditions_<kind>.cpp``. Numerics: :doc:`/numerics/domain_boundary_conditions`.

Add a diagnostic
----------------

.. flow::
   :horizontal:

   pharein/diagnostics.py | quantity list, checks
   -
   DiagnosticsManager | type names per model
   -
   type writer | detail/types/<type>.hpp
   -
   pharesee reader | field_qties, run.py | /data_analysis/reading_data

* The quantity string, e.g. ``/ions/pop/protons/density``, links everything; it becomes the file name
  ``ions_pop_protons_density.h5`` (:src:`src/diagnostic/detail/h5writer.hpp#L105-L111`).

New quantity in an existing type (e.g. ``fluid``):

#. Python: add it to ``FluidDiagnostics_.fluid_quantities``; check the ``population_name`` rules.
#. C++ (HDF5): list it in ``createFiles`` for the right tree, add branches in ``getDataSetInfo``,
   ``initDataSets``, ``write`` (``compute`` if needed) in ``detail/types/fluid.hpp``. Missing from
   ``createFiles``: ``Unknown Diagnostic Quantity`` (:src:`src/diagnostic/detail/h5typewriter.hpp#L143-L149`).
#. VTKHDF (optional): branches in ``setup`` / ``write`` of ``detail/vtk_types/fluid.hpp``.
#. Reader: ``field_qties`` (``pharesee/hierarchy/hierarchy_utils.py``), ``quantities_per_file`` and a
   ``Get...`` in ``pharesee/run/run.py``.
#. Test: add it to ``dump_all_diags`` in ``tests/diagnostic/__init__.py``.

New type ``foo``: a ``FooDiagnostics`` class naming blocks ``"foo" + count`` (C++ reads ``foo0``,
``foo1``..., :src:`src/diagnostic/diagnostic_manager.hpp#L51-L58`); the type name in ``registerDiagnostics``
(:src:`src/diagnostic/diagnostic_manager.hpp#L31-L44`), otherwise silently ignored; a writer in
``detail/types/foo.hpp``; a ``typeWriters_`` entry in ``h5writer.hpp`` *and* ``vtkh5_writer.hpp`` for each
model (``NullTypeWriter`` if unsupported).

Add a hybrid field
------------------

#. **Name:** an entry in ``HybridQuantity::Scalar`` (and ``Vector``, ``vector_components`` for a vector)
   (:src:`src/core/models/quantities/hybrid_quantities.hpp#L15-L50`).
#. **Centre:** its primal/dual centering in ``initLayoutCentering_``
   (:src:`src/core/data/grid/impl/yee/gridlayout_hybrid_yee.hpp#L52-L88`), and in ``yee_centering``
   (``pyphare/pyphare/core/gridlayout.py``) if pharesee reads it. See :doc:`gridlayout`.
#. **Own:** a ``Field`` / ``VecField`` member with a unique name, returned by
   ``getCompileTimeResourcesViewList`` or registered and allocated by its solver
   (``SolverPPC::registerResources``). See :doc:`resources_manager`.
#. **Use:** ``auto _ = resourcesManager->setOnPatch(*patch, myField);``.
#. **Ghosts:** a list in ``HybridMessengerInfo``, filled in ``fillMessengerInfo``, a refiner in
   ``registerQuantities``, and the fill called by the solver (e.g. ``fillFluxBorders``,
   :src:`src/amr/solvers/solver_ppc.hpp#L620-L623`). See :doc:`messengers`.
#. **Restarts:** nothing to do, every registered resource is saved
   (:src:`src/amr/resources_manager/resources_manager.hpp#L349-L355`).

**Pitfall.** Nothing checks that names are unique: a second resource with the same name is skipped and
both share the same data (``registerResource_``,
:src:`src/amr/resources_manager/resources_manager.hpp#L571-L587`).

Go deeper
---------

* :doc:`resources_manager` · :doc:`messengers` · :doc:`gridlayout` · :doc:`diagnostics_restarts`

Related
-------

* :doc:`working_on_phare` — build and tests.
