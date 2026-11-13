
========================
Adaptive Mesh Refinement
========================



Patch based approach
--------------------



Recursive time integration
--------------------------


Regridding
----------

Regridding is the process of changing the geometry of an existing mesh level
in the AMR hierarchy, by changing the number and geometry of some patches of
that level. This typically occurs when the needs for refinement evolve with
the solution over time.

Concretely, regridding means removing an entire level from the AMR hierarchy
and replacing it by a new one. The so-called "new" level is initialized by
copying values from the old level wherever they overlap, and by refining
values from coarser levels where they do not.

Technically, regridding involves different methods depending on whether
particle data or field data is refined. The specific operation performed also
depends on the type of field being refined: for instance, the magnetic field
refinement must preserve the divergence-free character of the field.


Field refinement
----------------


Particle refinement
-------------------



Field coarsening
----------------

Coarsening, also sometimes called "restriction", is the process of projecting
the solution existing on a given level onto the region it overlaps on the next
coarser level. Coarsening is used so that the solution on the coarse level
"feels" the fine solution in the overlapped region, which is assumed to be of
better quality.


Fields at level boundaries
--------------------------


Particle at level boundaries
----------------------------


.. _tagging-for-refinement:

Tagging for refinement
-----------------------


.. _clustering:

Clustering
----------


.. _tiling:

Tiling
~~~~~~


.. _berger-rigoutsos:

Berger-Rigoutsos
~~~~~~~~~~~~~~~~



