.. _dev-resources-manager:

================================
Resources manager and patch data
================================

.. container:: breadcrumb

   You are here: :doc:`/dev/index` › Resources manager

**In one sentence:** the *resources manager* is the adapter that registers the arrays
of PHARE's physics objects with :term:`SAMRAI`, and later points those objects at the
array of one :term:`patch`, so that the physics code never sees SAMRAI.

Why it exists
-------------

``Field``, ``VecField``, ``Ions``, ``Electromag``... live in ``src/core`` and know
nothing about SAMRAI. SAMRAI stores data per patch, as objects called *patch data*, and
patches come and go when the mesh is refined. So the model state is built once, as a set
of empty *views*, and the resources manager points them at the right patch data when
needed (``ResourcesManager``,
:src:`src/amr/resources_manager/resources_manager.hpp#L71-L113`).

Only two kinds of data are registered: ``Field`` (and the ``TensorField`` made of
fields) and ``ParticleArray`` packs
(:src:`src/amr/resources_manager/resources_manager.hpp#L80-L84`).

The picture
-----------

.. mermaid::

   flowchart TD
     C["core object<br/>Ions, Electromag, VecField..."] -->|"name, getCompileTimeResourcesViewList"| RM["ResourcesManager<br/>name to variable and id"]
     RM -->|"patch.getPatchData(variable, context)"| PD["SAMRAI patch data<br/>FieldData, ParticlesData"]
     PD -->|"getPointer()"| RM
     RM -->|"setBuffer(pointer)"| C
     click C href "particles_ions.html" "Particles and ions"
     click PD href "fields.html" "Field, VecField, tensor fields"

How it works
------------

#. **Discovery.** A core object says what it holds through two optional functions.
   ``getCompileTimeResourcesViewList()`` returns a tuple of references, for example
   ``(E, B)`` for ``Electromag``
   (:src:`src/core/data/electromag/electromag.hpp#L59-L64`);
   ``getRunTimeResourcesViewList()`` returns a vector, for example the populations of
   ``Ions`` (:src:`src/core/data/ions/ions.hpp#L261-L271`). A leaf is recognised by
   traits: ``is_field`` (has ``physicalQuantity()``), ``is_tensor_field`` (has
   ``components()``), ``is_particles`` (has ``setBuffer`` taking a ``ParticlesPack``)
   (:src:`src/amr/resources_manager/resources_manager_utilities.hpp#L23-L85`).
   Every function below walks this tree recursively until it reaches leaves
   (``registerResources``, :src:`src/amr/resources_manager/resources_manager.hpp#L169-L199`).

#. **Registration, once.** For each leaf, ``registerResource_`` builds a SAMRAI
   *variable* (what the data is) and registers it with a SAMRAI *context*
   (a label under which variables are looked up; one context, ``"default"``). SAMRAI
   returns an integer id. Name, variable and id are stored in a map keyed by the
   leaf's name (``registerResource_``,
   :src:`src/amr/resources_manager/resources_manager.hpp#L571-L588`).
   Registering a name already in the map does nothing. Solvers register their own
   temporaries the same way, for example ``SolverPPC::registerResources``
   (:src:`src/amr/solvers/solver_ppc.hpp#L211-L215`); the integrator calls it when a
   solver is registered
   (:src:`src/amr/multiphysics_integrator.hpp#L213`).

#. **Allocation, per patch.** When a patch is created, the integrator calls
   ``allocate`` of the model, the solver and the messenger
   (:src:`src/amr/multiphysics_integrator.hpp#L316-L318`).
   ``allocate_`` asks SAMRAI to create the patch data if the patch does not have it yet
   (:src:`src/amr/resources_manager/resources_manager.hpp#L626-L646`).

#. **Pointing a view at a patch.** ``setOnPatch(patch, obj1, obj2...)`` looks each
   leaf up by name, gets its patch data from the patch, and calls
   ``obj.setBuffer(pointer)`` (``setResourcesInternal_``,
   :src:`src/amr/resources_manager/resources_manager.hpp#L595-L609`). It returns a
   ``ResourcesGuard`` that undoes this when it goes out of scope
   (:src:`src/amr/resources_manager/resources_guards.hpp#L27-L58`).

#. **Looping on patches.** ``enumerate(level, obj...)`` returns a range over the
   patches of a level. Moving to a patch sets all listed objects on it; leaving it
   unsets them (``LevelLooper``,
   :src:`src/amr/resources_manager/resources_manager.hpp#L383-L454`). Prefer it to a
   hand-written loop. It was added in #1014 (resource manager ``enumerate``).

An example, shortened from a debug helper of the solver
(``_debug_log_move_ions``, :src:`src/amr/solvers/solver_ppc.hpp#L565-L579`):

.. code-block:: cpp

   for (auto& patch : rm.enumerate(level, model.state.ions))
   {
       for (auto& pop : model.state.ions)
           nbrLevelGhostParticles += pop.levelGhostParticles().size();
   }

Inside the loop body, ``model.state.ions`` and all its fields and particle arrays
point at the data of ``*patch``.

What the SAMRAI side holds
--------------------------

* ``FieldData`` derives from SAMRAI's ``PatchData`` and holds both a ``GridLayout``
  and the ``Grid`` (the memory) of the field
  (:src:`src/amr/data/field/field_data.hpp#L25-L52`,
  :src:`src/amr/data/field/field_data.hpp#L337-L338`). ``getPointer()`` returns the
  address of that ``Grid`` (:src:`src/amr/data/field/field_data.hpp#L294`), which is
  what ``Field::setBuffer`` reads
  (:src:`src/core/data/field/field.hpp#L59-L68`). The :term:`GridLayout` comes from
  the patch geometry (``FieldDataFactory::allocate``,
  :src:`src/amr/data/field/field_data_factory.hpp#L67-L80`).
* ``FieldVariable`` creates the factory that allocates a ``FieldData``
  (:src:`src/amr/data/field/field_variable.hpp#L38-L48`).
* ``ParticlesData`` holds five particle arrays (domain, patch ghost, level ghost, old and
  new level ghost) and a ``ParticlesPack`` of pointers to them
  (:src:`src/amr/data/particles/particles_data.hpp#L460-L468`); ``getPointer()``
  returns that pack (:src:`src/amr/data/particles/particles_data.hpp#L453`).
* ``TensorFieldData`` holds the grids of all components and returns them from
  ``getPointer()`` (:src:`src/amr/data/tensorfield/tensor_field_data.hpp#L319`);
  ``TensorField::setBuffer`` forwards component ``i`` to component ``i``
  (:src:`src/core/data/tensorfield/tensorfield.hpp#L84-L96`).
* The patch data types are chosen at compile time from the leaf's kind
  (``UserFieldType``, ``UserParticleType``, ``UserTensorFieldType``,
  :src:`src/amr/resources_manager/field_resource.hpp#L15-L19`,
  :src:`src/amr/resources_manager/particle_resource.hpp#L16-L24`,
  :src:`src/amr/resources_manager/tensor_field_resource.hpp#L17-L21`).

Times and restarts
------------------

SAMRAI keeps a time stamp on each patch data and uses it when it interpolates in time
between two stored states. ``setTime(obj, patch, time)`` sets it for every leaf of
``obj``; ``getTimes`` reads it
(:src:`src/amr/resources_manager/resources_manager.hpp#L260-L295`). For example the
messenger stamps the saved old fields ``Bold_`` and ``Eold_`` with the current time
(:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L591-L592`).

For :term:`restart`, ``registerForRestarts`` hands the ids of every registered resource
to SAMRAI's restart manager (:src:`src/amr/resources_manager/resources_manager.hpp#L328-L333`,
``restart_patch_data_ids`` L349-L356). ``FieldData`` and ``ParticlesData`` write and
read their arrays in ``putToRestart`` / ``getFromRestart``
(:src:`src/amr/data/field/field_data.hpp#L86-L97`,
:src:`src/amr/data/particles/particles_data.hpp#L144-L186`).

The tricks
----------

A guard unsets to ``nullptr``, not to the previous state
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** Two guards on the same object, one inside the other, are not safe.

**What the code does.** On exit, ``ResourcesGuard`` calls ``unsetResources_``, which
calls ``setBuffer(nullptr)`` on each object
(:src:`src/amr/resources_manager/resources_guards.hpp#L53-L58`,
``unsetResourcesInternal_``,
:src:`src/amr/resources_manager/resources_manager.hpp#L611-L622`). It does not restore
what was set before; the comment says so and lists "the active mechanism" as a TODO
(:src:`src/amr/resources_manager/resources_guards.hpp#L18-L25`). So if an inner scope
sets an object that an outer scope already set, the object is empty again when the inner
scope ends. There is no lock or check against nesting in the resources manager.

A forgotten guard is caught only in debug builds
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** The result of ``setOnPatch`` is easy to drop: a guard that is not kept is
destroyed at once and the views are empty again.

**Fix.** ``setOnPatch`` is marked ``NO_DISCARD``, so the compiler warns
(:src:`src/amr/resources_manager/resources_manager.hpp#L250-L256`). ``Field::operator()``
throws "Field is not usable" when the buffer is null, but only inside
``PHARE_DEBUG_DO``, that is in debug builds
(:src:`src/core/data/field/field.hpp#L74-L81`).

Resources are found by name
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** The map key is the object's name, so two objects with the same name share
one resource, and an unregistered name stops the run.

**Fix.** Registering an existing name is ignored
(:src:`src/amr/resources_manager/resources_manager.hpp#L579`); an empty or unknown name
throws (:src:`src/amr/resources_manager/resources_manager.hpp#L601-L606`).

.. admonition:: Known issue
   :class: note

   Restarts register every resource, including temporaries that need not be saved:
   see the comment pointing to issue #664 "Restarts writing redundant information"
   (:src:`src/amr/resources_manager/resources_manager.hpp#L349-L356`).

Go deeper
---------

* :doc:`/dev/fields`
* :doc:`/dev/gridlayout`

Related
-------

* :doc:`/dev/messengers`
* :doc:`/dev/particles_ions`
