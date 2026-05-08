
========================
Adaptive Mesh Refinement
========================

Introduction
------------

Hybrid PIC simulations of space plasmas must simultaneously capture large-scale dynamics and
small-scale kinetic structures such as current sheets, magnetic reconnection sites, and shock
transition layers. Resolving these structures everywhere in the domain is prohibitively expensive:
halving the cell size in 3D multiplies the number of cells by eight and, because the time step
scales with cell size, the total cost grows by a factor of sixteen.

Adaptive Mesh Refinement (AMR) addresses this by concentrating resolution only where it is needed.
A hierarchy of nested grids is maintained at runtime; fine grids cover the regions of interest while
the rest of the domain stays at coarse resolution. As structures move or develop, the fine grids
follow them.

Patch based approach
--------------------

PHARE uses the `SAMRAI <https://computing.llnl.gov/projects/samrai>`_ library as its AMR backbone.
The simulation domain is divided into rectangular **patches**, each storing a contiguous block of
field and particle data. Patches are grouped into **levels**:

- **Level 0** (coarsest) covers the entire domain and always exists.
- **Level 1, 2, …** cover sub-regions at progressively finer resolution.
- The **refinement ratio** between consecutive levels is fixed at **2**: the cell size halves at
  each level. This value is not user-configurable.
- Within a level, patches are non-overlapping; a cell belongs to exactly one patch.

The resulting structure is a patch hierarchy. SAMRAI handles the bookkeeping of which patches exist,
how they overlap across levels, and how data moves between them.

Two refinement strategies
-------------------------

PHARE supports two ways of deciding which regions to refine.

**Static boxes** (``refinement="boxes"``)
    The user provides a dictionary of boxes at configuration time. The hierarchy is fixed for the
    entire simulation. This is the simplest option and is useful when the regions of interest are
    known in advance.

    .. code-block:: python

        refinement_boxes = {"L0": {"B0": [(10,), (50,)]}}

    Each entry maps a level label (``"L0"`` = level 0) to named boxes defined by their lower and
    upper corner indices. The fine level covers exactly these boxes throughout the run.

**Adaptive tagging** (``refinement="tagging"``)
    At runtime, PHARE evaluates the magnitude of the magnetic field gradient in every coarse cell.
    Cells whose gradient exceeds ``tagging_threshold`` (default ``0.1``) are tagged, and SAMRAI
    creates or removes fine patches to cover the tagged region. As the plasma evolves and structures
    move, the fine grid tracks them automatically.

Recursive time integration
--------------------------

SAMRAI drives a recursive subcycling algorithm. Because the CFL condition ties the time step to the
cell size (and electromagnetic wave speeds scale with ``1/dx``), a level that is twice as fine must
take twice as many steps to advance by the same physical time. For diffusive or wave-dominated
problems the constraint is even tighter and the ratio is the square of the refinement ratio; with
refinement ratio 2 each fine level takes **4 substeps** for every coarse step
(see :doc:`temporal_discretization`).

The recursion proceeds as follows for a two-level hierarchy:

1. Level 0 advances one step ``dt₀``.
2. Level 1 advances four substeps of ``dt₁ = dt₀ / 4``, synchronizing ghost data with level 0 at
   each substep.
3. After level 1 has caught up to level 0, fine-to-coarse coarsening is applied (see
   `Field coarsening`_ below).
4. The process repeats for deeper levels before returning control to the coarser level.

Field refinement
----------------

When a new fine level is created, or when ghost cells on an existing fine level need to be filled,
field values are **interpolated from the coarser level**. PHARE uses spatial interpolation whose
order is set by the global ``interp_order`` parameter.

During subcycling, the coarse level has field snapshots at the beginning and end of the coarse step.
Ghost cells on the fine level at an intermediate time are filled by **time interpolation** between
these two snapshots, ensuring that the fine-level evolution sees a consistent coarse-level boundary
condition at every substep.

Particle refinement
-------------------

Particles do not live on the coarser level inside a refined region; they exist only on the finest
covering level. When a fine level is created over a coarse region, coarse particles that fall inside
are **split** into multiple fine particles.

The number of fine particles produced per coarse particle is controlled by
``refined_particle_nbr``. The splitting conserves moments: particle weights are adjusted so that
the zeroth moment (number density) and the first moment (momentum density) are preserved by the
collection of fine particles.

At level boundaries, a layer of **ghost particles** is populated from the coarser level using the
same splitting procedure, providing the fine level with a consistent particle environment just
outside its domain (see `Particles at level boundaries`_ below).

Field coarsening
----------------

After the fine level has advanced and is synchronized with the coarse level in time, fine field
data inside the refined region is **coarsened** back onto the coarse level. Each coarse cell
receives a weighted average of the fine cells it contains.

This step is essential for consistency: without it, the coarse-level solution inside a refined
region would lag behind the higher-resolution solution computed there. Coarsening ensures that
diagnostics and inter-level interactions always see the most accurate available data.

Fields at level boundaries
--------------------------

Ghost cells at the boundaries of a patch or level are filled differently depending on the source:

- **Same-level neighbor patch**: a direct copy of the interior data from the adjacent patch. No
  interpolation is needed because both patches share the same cell size.
- **Coarser level**: spatial interpolation from the overlapping coarse cells, combined with time
  interpolation when the coarse and fine levels are not synchronized in time (subcycling case).

The ghost cell width is determined by the stencil requirements of the numerical operators and by
``interp_order``.

Particles at level boundaries
------------------------------

Two kinds of boundary particles appear at level edges:

- **Level ghost particles**: populated from the coarser level by particle splitting. They occupy
  the ghost region just outside the fine domain and participate in current deposition and force
  calculations, ensuring that the fine-level particles near the boundary see a correct environment.
- **Patch ghost particles**: exchanged between neighboring patches at the *same* level. This is a
  simple copy with no splitting; it mirrors the same-level field ghost fill.

After each push, particles that have crossed a patch or level boundary are migrated to the correct
owner, and ghost particle layers are repopulated.

Patch size and clustering
--------------------------

Several parameters control how SAMRAI clusters tagged cells into patches:

``smallest_patch_size``
    Minimum number of cells along each dimension per patch. The default depends on
    ``interp_order`` and is chosen so that ghost regions do not exceed the patch interior.

``largest_patch_size``
    Maximum number of cells along each dimension per patch. Keeping patches smaller improves
    load balance on many MPI ranks.

``nesting_buffer``
    Minimum gap in coarse cells between the boundary of a fine level and the boundary of the
    next coarser level. Default is ``0``. A positive value prevents fine patches from touching
    the domain boundary of their parent level.

``clustering``
    Algorithm used to group tagged cells into rectangular patches:

    - ``"tile"`` (default): divides the tagged region into a regular tile grid. Fast and
      predictable, but may include some untagged cells.
    - ``"berger"``: uses the Berger-Rigoutsos algorithm for a tighter fit around tagged cells.
      Produces fewer wasted cells at the cost of higher setup time.

For further details on how these parameters are set in a simulation script, see
:doc:`../usage/simulation_inputs`. For the time-step constraints imposed by subcycling, see
:doc:`temporal_discretization`.
