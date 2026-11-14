.. _dev-fields:

==============================
Field, VecField, tensor fields
==============================

.. container:: breadcrumb

   You are here: :doc:`/dev/index` › Field, VecField, tensor fields

**In one sentence:** a ``Field`` is a light *view* on one :term:`patch`'s array of one
physical quantity (name, quantity, pointer); the memory belongs to the patch data, and
the same view is pointed at each patch in turn.

The picture
-----------

.. plot::
   :include-source: False
   :caption: The model's view ``EM_B_x`` owns no memory. Inside a patch loop the resources
      manager points it, by name, at the array of the current patch.

   from phare_draw import *
   from matplotlib.patches import Rectangle, FancyBboxPatch

   fig, ax = strip(2, titles=["Loop on patch A: the view points to A's array",
                              "Next iteration, patch B: same view, B's array"])
   for k, a in enumerate(ax):
       for x0, name in [(0, "patch A"), (6, "patch B")]:
           patch(a, (x0, 0, x0 + 4, 1.4))
           text(a, x0 + 2, 2.1, name, color=C["patch_edge"], bold=True)
           for c in range(4):
               a.add_patch(Rectangle((x0 + 0.3 + c * 0.85, 0.45), 0.8, 0.6,
                                     fc="white", ec=C["grey"], lw=1, zorder=2))
           text(a, x0 + 2, 1.65, "Grid: owns the array", color=C["grey"])
       a.add_patch(FancyBboxPatch((3.4, -2.2), 3.2, 0.9, boxstyle="round,pad=0.02,rounding_size=0.15",
                                  fc=C["hi"], ec=C["B"], lw=1.8, zorder=2))
       text(a, 5.0, -1.75, "Field EM_B_x (view)", color=C["B"], bold=True)
       on, off = (2, 8) if k == 0 else (8, 2)
       a.update_datalim([(0, 2.45)])  # room for the patch labels under the title
       arrow(a, (5.0, -1.3), (off, 0.45), color=C["grey"], lw=1.2, ls="--", ms=10)
       arrow(a, (5.0, -1.3), (on, 0.45), color=C["B"], lw=2.4, ms=18)
   finish(fig, ax)

* ``Grid`` owns the array and can hand out a ``Field`` view of itself
  (:src:`src/core/data/grid/grid.hpp#L15-L26`); patch data (``FieldData``) holds a ``Grid``
  (:src:`src/amr/data/field/field_data.hpp#L338`).
* ``Field<dim, Quantity>`` is an ``NdArrayView`` plus a name and a quantity
  (:src:`src/core/data/field/field.hpp#L24-L100`); copying it copies the pointer.
* ``resourcesManager->enumerate(level, B, E, ...)`` sets every listed view for each patch;
  in a range-for loop they are unset when the loop ends
  (:src:`src/amr/resources_manager/resources_manager.hpp#L373-L431`).
* ``VecField`` = ``TensorField`` of rank 1 (x, y, z); rank 2 is a symmetric tensor, six
  components ``_xx … _zz`` (:src:`src/core/data/tensorfield/tensorfield.hpp#L16-L35`). Each
  component has its own quantity, hence its own centering.

The tricks
----------

The name is the key
   **Problem.** Several objects read the same data (the ion density: electrons, Ohm,
   diagnostics).
   **Fix.** Views are matched to patch data by name: same name, same array
   (``setResourcesInternal_``, :src:`src/amr/resources_manager/resources_manager.hpp#L595-L609`).

Usable or settable
   **Problem.** A view used outside its patch loop reads a null pointer.
   **Fix.** ``isUsable()`` (has data) / ``isSettable()`` (no data)
   (:src:`src/core/data/field/field.hpp#L70-L71`); debug builds throw
   ``Field is not usable: <name>`` on access (:src:`src/core/data/field/field.hpp#L73-L81`).

Derivatives are views too
   **Problem.** Hand-written ``curl`` was copied in several operators.
   **Fix.** ``partial``, ``curl``, ``grad``, ``div`` return views that compute the derivative
   where read, at the right centering (:src:`src/core/numerics/differential_operators.hpp`).
   :pr:`1468` (shared curl, grad, div).

One scratch buffer
   **Problem.** Temporary arrays per class waste memory.
   **Fix.** One shared buffer per linked ``phare_core``; ``ScratchWrite`` throws if anyone
   else holds it, ``ScratchRead`` does not check for a writer
   (:src:`src/core/data/field/field_scratch.hpp#L11-L55`).

Go deeper
---------

* :doc:`/dev/gridlayout` — where each component lives, how its array is indexed.
* :doc:`/dev/resources_manager` — registration, allocation, ``enumerate``.

Related
-------

* :doc:`/dev/electromagnetic` · :doc:`/dev/particles_ions`
