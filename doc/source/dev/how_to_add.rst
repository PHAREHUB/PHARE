.. _dev-how_to_add:

==========
How to add
==========

.. container:: breadcrumb

   You are here: :doc:`/dev/index` › :doc:`/dev/working_on_phare` ›
   How to add.

**In one sentence:** step-by-step recipes for the three changes people make
most often — a boundary condition, a diagnostic, a field — with every file
to touch and the test to copy.

Each recipe follows the data from the Python input to where C++ uses it.
Several of these steps, if skipped, do not fail at compile time but at run
time, with one of the error messages quoted below.

Add a boundary condition
------------------------

There are two levels. Pick the right one.

* A **boundary type** is what users write (``"open"``, ``"reflective"``,
  ``"super-magnetofast-inflow"``, ``"none"``). It is a *recipe*: for each
  quantity, which field condition fills the :term:`ghost cells <ghost cell>`. Most new
  boundaries are only a new recipe.
* A **field condition** is a building block that fills the ghost cells of
  one field: ``None``, ``Dirichlet``, ``Neumann``, ``Symmetric``,
  ``AntiSymmetric``, and two divergence-free variants for vectors. Add one
  only if none of these fits.

.. mermaid::

   flowchart TD
     U["Simulation(domain_boundaries={'xlower': {'type': 'open'}, ...})"] --> PB["pharein/boundary.py<br/>Boundary subclass, populate_dict"]
     PB --> D["dict: grid/domain_boundaries/&lt;loc&gt;/type<br/>(integer of BoundaryType)"]
     D --> BM["BoundaryManager"] --> BF["BoundaryFactory::create<br/>switch on BoundaryType"]
     BF --> FC["registerFieldCondition&lt;kind&gt;(quantity)"]
     FC --> USE["field refine patch strategy<br/>fills ghost cells at physical boundaries"]
     click USE href "ghost_cells.html" "Ghost cells"

.. note::

   At the pinned commit, physical (non-periodic) boundaries exist for the
   MHD model only. Python refuses them for a hybrid model
   (:src:`pyphare/pyphare/pharein/boundary.py#L187-L192`), the hybrid model
   builds an empty ``BoundaryManager``
   (:src:`src/amr/physical_models/hybrid_model.hpp#L94-L103`), and particles
   leaving through a physical boundary are not handled yet.

New boundary type
^^^^^^^^^^^^^^^^^

#. **C++ enum.** Add the name at the end of ``enum class BoundaryType``
   (:src:`src/core/boundary/boundary_defs.hpp#L16`).
#. **Python sees the enum.** Add ``.value("my_name", core::BoundaryType::MyName)``
   to the ``py::enum_<core::BoundaryType>`` binding
   (:src:`src/python3/cpp_etc.cpp#L165-L169`). The Python string is the
   user name with ``-`` replaced by ``_``. Without it, writing the dict
   fails with a ``ValueError`` listing the allowed names.
#. **Python class.** In ``pyphare/pyphare/pharein/boundary.py``, add a
   dataclass deriving from ``Boundary`` with a class attribute ``type``
   (the user string) and its parameters as fields; override
   ``populate_dict`` if it has parameters (model:
   ``SuperMagnetofastInflowBoundary``). Declare it *above*
   ``_type_to_class`` in ``boundary.py``: that dict is built from the
   classes defined before it
   (:src:`pyphare/pyphare/pharein/boundary.py#L150-L156`), so there is no
   list to edit, but a class defined below it or elsewhere is unknown.
#. **The recipe.** In ``BoundaryFactory::create``, add a ``case`` to the
   switch (:src:`src/core/boundary/boundary_factory.hpp#L76-L87`) and a
   ``register_my_name_conditions_`` function that, for *every* scalar and
   vector quantity it receives, calls
   ``boundary->template registerFieldCondition<FieldBoundaryConditionType::...>(quantity)``
   (models: ``register_open_conditions_``, ``register_inflow_conditions_``).
#. **Test.** Add a case to ``checkDispatch`` in
   ``tests/core/numerics/boundary_condition/test_field_boundary_conditions_boundary_factory.cpp``;
   for an end-to-end check copy
   ``tests/functional/mhd_shock/shock_with_boundaries.py``.

**Pitfall.** A quantity without a registered condition makes the run stop
with ``Field boundary condition not found.`` when ghost cells are filled.
This is why the ``none`` type registers a do-nothing condition for every
quantity (comment at :src:`src/core/boundary/boundary_factory.hpp#L100-L101`).

New field condition
^^^^^^^^^^^^^^^^^^^

#. Add an entry to ``FieldBoundaryConditionType``
   (``src/core/numerics/boundary_condition/field_boundary_condition.hpp``).
#. Write ``field_boundary_condition_<kind>.hpp`` next to it, deriving from
   ``IFieldBoundaryCondition`` (model: ``field_boundary_condition_neumann.hpp``).
#. Map the enum to the class with a ``FieldBoundaryConditionSelector``
   specialisation
   (:src:`src/core/numerics/boundary_condition/field_boundary_condition_resolver.hpp#L18-L75`).
#. Copy one of the ``test_field_boundary_conditions_<kind>.cpp`` tests and
   add the suffix to the list in
   ``tests/core/numerics/boundary_condition/CMakeLists.txt``.

See :doc:`/numerics/domain_boundary_conditions` for what each condition does
numerically.

Add a diagnostic
----------------

The *quantity string* (e.g. ``/ions/pop/protons/density``) is the single
key that links everything: the Python class sets it, the C++ type writer
tests it, it becomes the file name (``ions_pop_protons_density.h5``,
``H5Writer::fileString``, :src:`src/diagnostic/detail/h5writer.hpp#L105-L111`),
and the reader maps the file back.

.. mermaid::

   flowchart TD
     P["pharein/diagnostics.py<br/>quantity list + checks"] --> G["populate_diagnostics<br/>(generic, nothing to edit)"]
     G --> M["DiagnosticsManager<br/>registerDiagnostics: type list"]
     M --> W["type writer<br/>detail/types/&lt;type&gt;.hpp"]
     W --> F["&lt;quantity&gt;.h5"]
     F --> R["pharesee reader<br/>field_qties, quantities_per_file"]
     click M href "diagnostics_restarts.html" "Diagnostics internals"

New quantity in an existing type
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Example: a new ``fluid`` quantity.

#. **Python.** Add it to ``FluidDiagnostics_.fluid_quantities`` and check
   the rules in ``_setSubTypeAttributes`` (which quantities need or forbid
   a ``population_name``) in ``pyphare/pyphare/pharein/diagnostics.py``.
#. **C++ type writer.** In ``src/diagnostic/detail/types/fluid.hpp``, add the
   name to the list in ``createFiles`` for the right tree (``/ions/`` or
   ``/ions/pop/<pop>/``), and a branch in ``getDataSetInfo``,
   ``initDataSets`` and ``write`` (and ``compute`` if it must be computed
   before writing, like ``momentum_tensor``).
#. **Model view.** If the writer cannot reach the field yet, add an
   accessor in ``src/diagnostic/diagnostic_model_view.hpp``.
#. **VTKHDF (optional).** Same branches in
   ``src/diagnostic/detail/vtk_types/fluid.hpp``.
#. **Reader.** Map dataset names in ``field_qties``
   (``pyphare/pyphare/pharesee/hierarchy/hierarchy_utils.py``), and add the
   file to ``quantities_per_file`` with a ``Get...`` method in
   ``pyphare/pyphare/pharesee/run/run.py``.
#. **Test.** Add it to ``dump_all_diags`` in ``tests/diagnostic/__init__.py``,
   which the diagnostics and restart tests use.

**Pitfall.** If ``createFiles`` does not list the quantity, the run stops
with ``Unknown Diagnostic Quantity: <quantity>``
(``h5FileForQuantity``, :src:`src/diagnostic/detail/h5typewriter.hpp#L143-L149`).
Python and C++ lists are separate: check that both accept the same
combinations.

New diagnostic type
^^^^^^^^^^^^^^^^^^^

#. **Python.** A new ``FooDiagnostics(Diagnostics)`` class with its quantity
   list and ``type = "foo"`` in ``pharein/diagnostics.py``, exported in
   ``pharein/__init__.py``. Its ``__init__`` must name each block
   ``"foo" + count`` like the existing types: C++ looks for ``foo0``,
   ``foo1``... in order
   (:src:`src/diagnostic/diagnostic_manager.hpp#L51-L58`).
#. **Register the type name** in the model's list in ``registerDiagnostics``
   (:src:`src/diagnostic/diagnostic_manager.hpp#L31-L44`). A type not in
   this list is silently ignored.
#. **Writer.** A new ``src/diagnostic/detail/types/foo.hpp`` deriving from
   ``H5TypeWriter`` (models: ``fluid.hpp``, ``info.hpp``).
#. **Both formats.** Add ``{"foo", make_writer<...>()}`` to ``typeWriters_``
   in ``h5writer.hpp`` *and* in ``vtkh5_writer.hpp`` (``NullTypeWriter`` if
   VTKHDF does not support it), in the branch of each model (hybrid, MHD)
   that should offer the type. A missing entry throws when the type is
   looked up.
#. **Reader and tests** as above.

Add a hybrid field
------------------

A field that lives on every patch (like ``J`` or the electron pressure).

.. mermaid::

   flowchart TD
     Q["HybridQuantity entry"] --> C["centering table<br/>(C++ and pharesee)"]
     C --> O["Field member with a<br/>unique name"]
     O --> RM["registered in the<br/>resources manager"]
     RM --> U["setOnPatch: used<br/>in solver code"]
     RM --> G["ghosts: messenger<br/>refiner + fill"]
     RM --> RS["restarts: automatic"]
     click RM href "resources_manager.html" "Resources manager"
     click G href "messengers.html" "Messengers"

#. **Name it.** Add an entry to ``HybridQuantity::Scalar`` before ``count``;
   for a vector, also to ``Vector`` and to ``vector_components``
   (:src:`src/core/models/quantities/hybrid_quantities.hpp#L15-L50`).
#. **Centre it.** Give each component its primal/dual centering in
   ``initLayoutCentering_`` and add it to the table
   (:src:`src/core/data/grid/impl/yee/gridlayout_hybrid_yee.hpp#L52-L88`).
   If the field is written to disk and read back with ``pharesee``, add the
   same centering to ``yee_centering`` in
   ``pyphare/pyphare/core/gridlayout.py`` (one dict per direction). See
   :doc:`gridlayout`.
#. **Own it.** Declare a ``Field`` / ``VecField`` member with a *unique*
   name in the class that owns it (model state, solver, or model), and make
   it visible to the resources manager: return it from
   ``getCompileTimeResourcesViewList`` in a state class, or call
   ``registerResources`` and ``allocate`` explicitly in a solver
   (model: ``SolverPPC::registerResources`` and ``SolverPPC::allocate``).
   See :doc:`resources_manager`.
#. **Use it.** ``auto _ = resourcesManager->setOnPatch(*patch, myField);``
   makes it point to the patch's data.
#. **Ghosts (if needed).** Add a list to ``HybridMessengerInfo``, fill it in
   ``fillMessengerInfo``, build the refiner in
   ``HybridHybridMessengerStrategy::registerQuantities``, and call the
   matching fill from the solver at the right point of the step (e.g.
   ``fromCoarser.fillFluxBorders``,
   :src:`src/amr/solvers/solver_ppc.hpp#L620-L623`). See
   :doc:`messengers`.
#. **Restarts: nothing to do.** Every registered resource is saved
   (``restart_patch_data_ids``,
   :src:`src/amr/resources_manager/resources_manager.hpp#L349-L355`).
#. **Test.** Copy ``tests/amr/resources_manager/test_resources_manager.cpp``
   (register, allocate, ``setOnPatch``) and the gridlayout tests in
   ``tests/core/data/gridlayout/``.

**Pitfall.** Names must be unique, and nothing checks it: a second
resource with an already registered name is silently skipped, so both
objects end up pointing to the *same* data on the patch
(``ResourcesManager::registerResource_``,
:src:`src/amr/resources_manager/resources_manager.hpp#L571-L587`).

Go deeper
---------

* :doc:`resources_manager`, :doc:`messengers`, :doc:`gridlayout`
* :doc:`diagnostics_restarts`

Related
-------

* :doc:`working_on_phare` — build and tests.
