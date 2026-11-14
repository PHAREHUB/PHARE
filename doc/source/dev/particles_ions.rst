.. _dev-particles-ions:

================================
ParticleArray, Ions, populations
================================

.. container:: breadcrumb

   You are here: :doc:`/dev/index` › ParticleArray, Ions, populations

**In one sentence:** the ions of a :term:`hybrid model` are a list of :term:`populations <population>`
(species with their own mass and charge); each population holds its moments
(fields) and, on each :term:`patch`, five ``ParticleArray`` objects: its domain
particles and four kinds of ghost particles.

The picture
-----------

.. mermaid::

   flowchart TD
     I["<b>Ions</b><br/>total moments: charge density, mass density,<br/>bulk velocity, charge flux, momentum tensor"] --> P1["<b>IonPopulation</b> 'protons'<br/>mass, charge<br/>density, flux, momentum tensor"]
     I --> P2["<b>IonPopulation</b> 'alpha' ..."]
     P1 --> K["<b>ParticlesPack</b><br/>5 pointers, set per patch"]
     K --> D["domain particles"]
     K --> G["patch ghost particles"]
     K --> L["level ghost particles<br/>(pushable)"]
     K --> O["level ghost Old"]
     K --> N["level ghost New"]
     S["<b>ParticlesData</b><br/>SAMRAI patch data,<br/>owns the 5 ParticleArrays"] -. "pointers set by the<br/>resources manager" .-> K
     click G href "grabbing_particles.html" "Patch ghost particles"
     click L href "level_ghost_particles.html" "Level ghost particles"
     click O href "level_ghost_particles.html" "Level ghost particles"
     click N href "level_ghost_particles.html" "Level ghost particles"
     click S href "resources_manager.html" "Resources manager and patch data"

How it works
------------

One particle
^^^^^^^^^^^^

A ``Particle<dim>`` (:src:`src/core/data/particles/particle.hpp#L54-L60`) has:

* ``weight``: its density contribution (density per particle; see
  :doc:`/dev/particle_splitting` for how weights change with refinement);
* ``iCell``: the index of the cell it is in, in the level's global (:term:`AMR`) index
  space, so that it stays valid when the particle moves between patches;
* ``delta``: its position inside that cell, between 0 and 1 in each direction;
* ``v``: its three velocity components.

There is no charge or mass per particle: they belong to the population
(comment at L43; ``IonPopulation::charge``, ``mass``).

**History.** #1473 (particles carry no charge).

One population
^^^^^^^^^^^^^^

``IonPopulation``
(:src:`src/core/data/ions/ion_population/ion_population.hpp#L19-L146`) has a
name, a mass and a charge read from the input, its moments (particle density,
charge density, flux, momentum tensor) and a ``ParticlesPack``
(:src:`src/core/data/ions/ion_population/particle_pack.hpp#L19-L36`).

The pack holds five *pointers*, not arrays. The arrays themselves live in the
:term:`SAMRAI` patch data ``ParticlesData``
(:src:`src/amr/data/particles/particles_data.hpp#L460-L466`), one per patch.
Before working on a patch, the code asks the resources manager to point the
population's pack and fields at that patch's data (``setOnPatch``; see
:doc:`/dev/resources_manager`). The same ``IonPopulation`` object thus works on
every patch in turn.

.. list-table:: The five particle arrays of a patch
   :header-rows: 1
   :widths: 28 72

   * - Array
     - Holds
   * - ``domainParticles``
     - particles whose cell is in the patch and owned by it; the only ones
       that are really "the plasma" of this patch
   * - ``patchGhostParticles``
     - particles that just left the patch for a same-level neighbour; live
       between a push and the next exchange (:doc:`/dev/grabbing_particles`)
   * - ``levelGhostParticles``
     - the pushable copy of ``levelGhostParticlesOld``: fine particles split
       from the coarse level, just outside a refined level
       (:doc:`/dev/level_ghost_particles`)
   * - ``levelGhostParticlesOld``, ``levelGhostParticlesNew``
     - coarse particles split at the start and end of the coarse step, used for
       the moments near the level border

Ions
^^^^

``Ions`` (:src:`src/core/data/ions/ions.hpp#L295-L300`) holds the vector of
populations and the total moments, computed from the populations' moments:
charge density, mass density, bulk velocity, charge flux and the full momentum
tensor (``computeChargeDensity`` and the next functions, L85-L197; formulas on
:doc:`/dev/moments`).

ParticleArray: particles indexed by cell
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``ParticleArray<dim>``
(:src:`src/core/data/particles/particle_array.hpp#L21-L238`) is a
``std::vector`` of particles (an *array of structures*: each particle's fields
are stored together) plus a ``CellMap``: for every
cell of a box, the list of indexes of the particles in that cell
(:src:`src/core/utilities/cellmap.hpp#L26-L38`). The map is kept up to date
when particles are added, erased or change cell (``push_back``,
``erase`` on an index range, ``change_icell``, L99 and L131-L190). The other
``erase``, on a pair of iterators, does *not* update the map; the comment at
L101-L112 explains why this is considered harmless.

What it gives:

* ``nbr_particles_in(box)``: count the particles in a box by looking only at
  the cells of that box;
* ``export_particles(box, dest)``: copy the particles of a box into another
  array (used by the exchanges between patches);
* ``partition(range, predicate)``: move to the end of a range the particles
  whose *cell* fails a test, and return the range of those that pass
  (:src:`src/core/utilities/cellmap.hpp#L400-L445`). The pusher's *selectors*
  ("is it still in my domain?", "is it in the ghost layer?") are partitions.

The tricks
----------

Selections work on cells, not on particles
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** Testing every particle against a box is slow when only a few
cells matter, and the pusher sorts particles several times per step.

**Fix.** All selections go through the cell map: the test is done once per
cell, then applied to all the particles indexed in that cell. The catch: a
particle whose cell is outside the map's box is invisible to the selections.

The cell map covers more than the ghost box
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** A pushed particle can end up a few cells past the ghost box: the
pusher lets a particle's cell index change by up to 2 per half push (see
:doc:`/dev/grabbing_particles`). If its new cell were outside the map, the
selectors would not see it, and it would be kept or dropped wrongly.

**Fix.** Each array of ``ParticlesData`` maps the ghost box grown by a margin
(``ghostSafeMapLayer``,
:src:`src/amr/data/particles/particles_data.hpp#L68-L73` and L79-L87), equal to
``particleMapMargin`` = max(2, 4 − particle ghost width)
(:src:`src/core/data/grid/gridlayout.hpp#L59-L68`). The comment there gives the
reasoning: level ghosts are selected after each half push, domain particles
after both.

**History.** #1592 (selectors missed far particles).

One SAMRAI variable for five arrays
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** SAMRAI moves data per *variable*. Having five separate variables
per population would multiply the schedules.

**Fix.** ``ParticlesData`` is one SAMRAI *patch data* (the object SAMRAI stores per
patch for a variable) holding all five arrays. SAMRAI tells it what to copy
with an *overlap* object, built by a fill pattern (:doc:`/dev/fill_patterns`);
how each operation treats the arrays depends on the overlap's type:
a ``ParticlesDomainOverlap`` means "grab patch ghosts into my domain", a plain
cell overlap means "copy domain particles"
(:src:`src/amr/data/particles/particles_data.hpp#L245-L270`). Refinement goes
through the split operators instead (:doc:`/dev/particle_splitting`).

Restarts store four arrays out of five
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Restart files hold the domain particles and the three level ghost arrays, not
the patch ghost particles, which are empty between time steps
(``ParticlesData::putToRestart``,
:src:`src/amr/data/particles/particles_data.hpp#L114-L141`).

Go deeper
---------

* :doc:`/dev/grabbing_particles` — what happens to patch ghost particles.
* :doc:`/dev/level_ghost_particles` — the three level ghost arrays.
* :doc:`/dev/moments` — how particles become moments.

Related
-------

* :doc:`/dev/resources_manager` — how the pack is pointed at a patch's data.
* :doc:`/dev/ppc_order` — when particles are pushed and deposited.
* :doc:`/numerics/particles` — shapes, pusher, interpolation.
