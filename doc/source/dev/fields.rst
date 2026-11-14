.. _dev-fields:

==============================
Field, VecField, tensor fields
==============================

.. container:: breadcrumb

   You are here: :doc:`/dev/index` › Field, VecField, tensor fields

**In one sentence:** a ``Field`` is a light *view* on a patch's array of one physical
quantity: a name, the quantity (which fixes its centering), and a pointer to data that
someone else owns; ``VecField`` and ``TensorField`` group three or six of them.

Why it exists
-------------

The model state (``Electromag``, ``Ions``, the MHD state...) is built once for a whole
simulation, but the data live in :term:`SAMRAI`, one array per :term:`patch`, and patches come and go
when the mesh is refined. So the objects the solvers manipulate cannot own their data.
They are views that are pointed at one patch's data, used, then pointed at the next
patch. ``Field`` is that view for one scalar array.

The picture
-----------

.. mermaid::

   flowchart TD
     P["SAMRAI patch"] --> FD["FieldData<br/>owns a Grid (the memory)"]
     FD --> RM["ResourcesManager<br/>finds the data by name"]
     RM -->|"setBuffer(pointer)"| F["Field (view)<br/>name, quantity, pointer"]
     F --> VF["VecField / TensorField<br/>3 or 6 Fields"]
     VF --> S["model state<br/>B, E, J, ion moments"]
     click RM href "resources_manager.html" "Resources manager"
     click F href "fields.html#how-it-works" "How it works"
     click S href "electromagnetic.html" "Faraday, Ampère, Ohm"

How it works
------------

#. **Grid owns, Field views.** ``Grid`` owns an n-dimensional array and carries its name
   and quantity (:src:`src/core/data/grid/grid.hpp#L15-L26`). Its own comment states the
   intent: separate who allocates memory from who computes on it. A ``Grid`` hands out a
   ``Field`` view of itself (``operator&``, ``operator*``,
   :src:`src/core/data/grid/grid.hpp#L90-L94`). The SAMRAI patch data for a field,
   ``FieldData``, holds a ``Grid`` (:src:`src/amr/data/field/field_data.hpp#L338`).
#. **A Field is a view with a name.** ``Field<dim, Quantity>`` derives from
   ``NdArrayView``: a pointer and a shape, plus a name and a physical quantity
   (:src:`src/core/data/field/field.hpp#L24-L100`). Copying a ``Field`` copies the pointer,
   not the values.
#. **Pointing a view at a patch.** The resources manager looks up the patch data by the
   view's *name* and calls ``setBuffer`` with its pointer
   (``setResourcesInternal_``, :src:`src/amr/resources_manager/resources_manager.hpp#L595-L609`).
   Solvers do this through ``resourcesManager->enumerate(level, B, E, ...)``, which sets
   every listed view for each patch of the level in turn, and unsets them afterwards
   (details: :doc:`/dev/resources_manager`).
#. **Usable or settable.** A view with data is *usable*; a view without data is
   *settable* (``isUsable``, ``isSettable``,
   :src:`src/core/data/field/field.hpp#L70-L71`). Composite objects (``Ions``,
   ``Electrons``, ...) are usable only when all their views are. Operators check this
   before working.
#. **Vectors and tensors.** ``VecField`` is ``TensorField`` of rank 1: three ``Field``
   components x, y, z (:src:`src/core/data/vecfield/vecfield.hpp`). Rank 2 is a
   *symmetric* 3×3 tensor stored as six components named ``_xx, _xy, _xz, _yy, _yz, _zz``,
   used for the ion momentum tensor
   (:src:`src/core/data/tensorfield/tensorfield.hpp#L16-L35`). Each component has its own
   quantity, so its own centering (for B: ``Bx`` on x-faces, ``By`` on y-faces, ...).

The tricks
----------

The name is the key
^^^^^^^^^^^^^^^^^^^

**Problem.** Several objects need the same data: the ion density is read by the
electrons, Ohm's law and the diagnostics.

**Fix.** Views are matched to patch data by name. Two views with the same name, set on
the same patch, see the same array. ``Field::setBuffer`` checks the names agree in debug
builds (:src:`src/core/data/field/field.hpp#L59-L68`). A component of a ``VecField``
named ``EM_B`` is ``EM_B_x``, ``EM_B_y``, ``EM_B_z``
(:src:`src/core/data/tensorfield/tensorfield.hpp#L26-L35`).

Reading a view that has no data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** Using a view outside an ``enumerate`` loop reads through a null pointer.

**Fix.** In debug builds, ``Field::operator()`` throws ``Field is not usable: <name>``
(:src:`src/core/data/field/field.hpp#L73-L81`). In release builds there is no check, for
speed. If a debug run throws this, the view was used outside the patch loop that sets it.

Derivatives are views too
^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** ``curl B`` written out by hand, per dimension, was copied in several
operators.

**Fix.** ``partial``, ``curl``, ``grad``, ``div`` return small views that compute the
derivative when read, at the right centering (:src:`src/core/numerics/differential_operators.hpp`).
Nothing is allocated; ``Faraday`` reads ``curl(E)`` directly into its update.

**History.** #1468 (shared curl, grad, div).

One scratch buffer
^^^^^^^^^^^^^^^^^^

**Problem.** Some operators need a temporary array. Static buffers in each class
multiply memory and can be overwritten by nested calls.

**Fix.** There is one shared scratch buffer per linked ``phare_core`` library (one per process in practice). Code takes it through
``ScratchWrite`` (exclusive) or ``ScratchRead`` (shared, read-only); the holders are
counted so that a nested write cannot silently overwrite data still in use
(:src:`src/core/data/field/field_scratch.hpp#L11-L55`).

Go deeper
---------

* :doc:`/dev/gridlayout`: where each component lives and how its arrays are indexed.
* :doc:`/dev/resources_manager`: registration, allocation and ``enumerate``.

Related
-------

* :doc:`/dev/electromagnetic`: the operators that use these fields.
* :doc:`/dev/particles_ions`: ions and populations, built from the same views.
