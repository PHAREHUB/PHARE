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

* Core objects know nothing about SAMRAI, which stores data per patch (*patch data*).
* The model state is built once, as empty *views*; the resources manager points them
  at the right patch data
  (:src:`src/amr/resources_manager/resources_manager.hpp#L71-L113`).
* Registered: ``Field``, ``TensorField``, ``ParticleArray`` packs
  (:src:`src/amr/resources_manager/resources_manager.hpp#L80-L84`).

The picture
-----------

.. plot::
   :include-source: false
   :alt: Registration once, then for each patch the resources manager points the views of the core object at the patch data.

   from phare_draw import *
   from matplotlib.patches import FancyBboxPatch

   fig, (a, b) = strip(2, ["registration, once", "setOnPatch(patch 1), for each patch"], height=2.4)

   def box(ax, x, y, w, h, s, ec=None, fc=None):
       ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.12",
                                   fc=fc or C["faint"], ec=ec or C["grey"], lw=1.6, zorder=2))
       ax.text(x + w / 2, y + h / 2, s, ha="center", va="center", fontsize=11, zorder=3)

   X = (0, 4.8, 9.6)   # box left edges, width 2.6
   box(a, X[0], 0, 2.6, 1.0, "core object\nviews")
   box(a, X[1], 0, 2.6, 1.0, "Resources-\nManager")
   box(a, X[2], 0, 2.6, 1.0, "SAMRAI\nvariable, id")
   arrow(a, (2.65, 0.5), (4.75, 0.5), color=C["grey"], lw=1.6, text="name", text_off=(0, 0.3))
   arrow(a, (7.45, 0.5), (9.55, 0.5), color=C["grey"], lw=1.6, text="register", text_off=(0, 0.3))
   a.plot([0, 12.2], [-0.4, 1.4], alpha=0)

   box(b, X[0], 0.9, 2.6, 1.0, "core object\nviews (empty)")
   box(b, X[1], 0.9, 2.6, 1.0, "Resources-\nManager", ec=C["L1"])
   box(b, X[2], 1.5, 2.6, 0.8, "patch 0 data")
   box(b, X[2], 0.5, 2.6, 0.8, "patch 1 data", ec=C["L1"], fc=C["patch"])
   arrow(b, (7.45, 1.2), (9.55, 0.9), color=C["L1"], lw=1.8, text="getPointer()", text_off=(-0.1, -0.45))
   arrow(b, (4.75, 1.4), (2.65, 1.4), color=C["L1"], lw=1.8, text="setBuffer(ptr)", text_off=(0, 0.3))
   b.plot([0, 12.2], [0.1, 2.5], alpha=0)
   finish(fig, [a, b])

Look at the bottom row: the pointer of one patch's data goes to the views.

How it works
------------

* **Discovery.** An object lists its parts with ``getCompileTimeResourcesViewList()``
  (:src:`src/core/data/electromag/electromag.hpp#L59-L64`) or
  ``getRunTimeResourcesViewList()``
  (:src:`src/core/data/ions/ions.hpp#L261-L271`); ``registerResources`` walks the tree
  (:src:`src/amr/resources_manager/resources_manager.hpp#L169-L199`).
* **Registration, once.** ``registerResource_`` builds a SAMRAI variable; the id is kept by name
  (:src:`src/amr/resources_manager/resources_manager.hpp#L571-L588`).
* **Allocation.** ``allocate_`` creates missing patch data
  (:src:`src/amr/resources_manager/resources_manager.hpp#L626-L646`).
* **Pointing.** ``setOnPatch(patch, obj...)`` calls ``obj.setBuffer(pointer)`` and
  returns a ``ResourcesGuard``
  (:src:`src/amr/resources_manager/resources_manager.hpp#L595-L609`,
  :src:`src/amr/resources_manager/resources_guards.hpp#L27-L58`).
* **Loops.** Prefer ``enumerate(level, obj...)``
  (:src:`src/amr/resources_manager/resources_manager.hpp#L383-L454`; :pr:`1014`).

What the SAMRAI side holds
--------------------------

* ``FieldData``: a ``GridLayout`` and the ``Grid``; ``getPointer()`` returns it
  (:src:`src/amr/data/field/field_data.hpp#L25-L52`,
  :src:`src/amr/data/field/field_data.hpp#L294`).
* ``ParticlesData``: five particle arrays and a ``ParticlesPack`` of pointers to them
  (:src:`src/amr/data/particles/particles_data.hpp#L460-L468`).
* ``TensorFieldData``: the grids of all components
  (:src:`src/amr/data/tensorfield/tensor_field_data.hpp#L319`).

Times and restarts
------------------

* Each patch data carries a time
  (:src:`src/amr/resources_manager/resources_manager.hpp#L260-L295`).
* :term:`restart`: ``registerForRestarts`` hands every id to SAMRAI
  (:src:`src/amr/resources_manager/resources_manager.hpp#L328-L333`).

The tricks
----------

A guard unsets to ``nullptr``, not to the previous state
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** Nested guards on the same object are not safe.

**What the code does.** On exit, ``ResourcesGuard`` calls ``setBuffer(nullptr)``,
not the previous state; nothing checks nesting
(:src:`src/amr/resources_manager/resources_guards.hpp#L18-L25`,
:src:`src/amr/resources_manager/resources_manager.hpp#L611-L622`).

A forgotten guard is caught only in debug builds
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** A guard not kept is destroyed at once; the views are empty.

**Fix.** ``setOnPatch`` is ``NO_DISCARD``
(:src:`src/amr/resources_manager/resources_manager.hpp#L250-L256`);
``Field::operator()`` throws only in debug builds
(:src:`src/core/data/field/field.hpp#L74-L81`).

Resources are found by name
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** Same name, same resource; an unknown name stops the run.

**Fix.** Re-registering a name is ignored
(:src:`src/amr/resources_manager/resources_manager.hpp#L579`); an empty or unknown name
throws (:src:`src/amr/resources_manager/resources_manager.hpp#L601-L606`).

.. admonition:: Known issue
   :class: note

   Restarts register every resource, including temporaries (issue :issue:`664`,
   :src:`src/amr/resources_manager/resources_manager.hpp#L349-L356`).

Go deeper
---------

* :doc:`/dev/fields`
* :doc:`/dev/gridlayout`

Related
-------

* :doc:`/dev/messengers`
* :doc:`/dev/particles_ions`
