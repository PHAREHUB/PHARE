.. _dev-particles-ions:

================================
ParticleArray, Ions, populations
================================

.. container:: breadcrumb

   You are here: :doc:`/dev/index` › ParticleArray, Ions, populations

**In one sentence:** the ions of a :term:`hybrid model` are a list of
:term:`populations <population>` (species with their own mass and charge); each
population holds its moments and, on each :term:`patch`, five ``ParticleArray``
objects: its domain particles and four kinds of ghost particles.

* ``Ions`` → vector of ``IonPopulation`` + total moments (charge and mass
  density, bulk velocity, charge flux, momentum tensor)
  (:src:`src/core/data/ions/ions.hpp#L295-L300`).
* ``IonPopulation`` → name, mass, charge, density, flux, momentum tensor, and a
  ``ParticlesPack`` of 5 *pointers*
  (:src:`src/core/data/ions/ion_population/ion_population.hpp#L137-L144`).
* ``ParticlesData`` (SAMRAI patch data, one per patch) owns the 5 arrays; the
  resources manager points the pack at them on each patch
  (:src:`src/amr/data/particles/particles_data.hpp#L460-L466`,
  :doc:`/dev/resources_manager`).

The five arrays of a patch
--------------------------

.. plot::
   :include-source: False
   :caption: Where the particles of each array sit, for patch A on a refined
      level whose neighbour B is on the right (2 ghost cells).

   from phare_draw import *
   from matplotlib.patches import Rectangle
   fig, (a,) = strip(1, ["Patch A and its five particle arrays"], height=3.6)
   grid(a, (-2, -2, 13, 8), color=C["light"])
   patch(a, (0, 0, 8, 6), ghosts=2, label="patch A")
   patch(a, (8, 0, 13, 6), label="patch B")
   a.add_patch(Rectangle((8, 0), 2, 6, fc=C["E"], alpha=0.3, ec="none", zorder=1))
   particles(a, [(1.5, 1.4), (3.2, 4.1), (5.6, 2.3), (6.8, 4.6), (2.4, 2.8)], color=C["ion"])
   particles(a, [(8.6, 3.3), (9.4, 1.2)], color=C["L2"])
   particles(a, [(-1.4, 2.2), (2.3, 6.8), (5.1, -1.3)], color=C["grey"])
   particles(a, [(-0.7, 3.6), (4.2, 7.3), (6.5, -0.6)], color=C["L1"])
   text(a, 4, 3.0, "domain", color=C["ion"], bold=True)
   text(a, 11.2, 3.6, "patch ghost", color=C["L2"], bold=True)
   text(a, 4.0, -2.8, "level ghost: Old (grey), New (blue), copy of Old", color=C["ghost_edge"], bold=True)
   a.plot([4.0], [-3.4], alpha=0)
   finish(fig, [a])

* ``domainParticles``: cell in the patch and owned by it — the plasma of A.
* ``patchGhostParticles``: just left A for B; deposited by A, then grabbed by B
  (:doc:`/dev/grabbing_particles`).
* ``levelGhostParticlesOld`` / ``New`` / ``levelGhostParticles`` (pushable copy):
  split from the coarse level, outside the level only (:doc:`/dev/level_ghost_particles`).
* A ``Particle``: ``weight`` (density per particle), ``iCell`` (global AMR
  cell index), ``delta`` (position in the cell, in [0,1)), ``v``; no charge
  (:src:`src/core/data/particles/particle.hpp#L54-L60`). :pr:`1473`.

ParticleArray: a vector and a cell map
--------------------------------------

.. plot::
   :include-source: False
   :caption: A ParticleArray is a vector of particles plus, for each cell, the
      list of indexes of its particles. ``partition`` tests cells, not
      particles, and moves the particles of failing cells to the end.

   from phare_draw import *
   from matplotlib.patches import Rectangle
   cells = {0: "c0", 1: "c1", 2: "c2", 3: "c3"}
   col = {0: C["L1"], 1: C["ion"], 2: C["rho"], 3: C["bad"]}
   before = [("p0", 1), ("p1", 3), ("p2", 0), ("p3", 2), ("p4", 3), ("p5", 1)]
   after = [("p0", 1), ("p5", 1), ("p2", 0), ("p3", 2), ("p4", 3), ("p1", 3)]
   def vec(a, items, y, title):
       text(a, -0.3, y + 0.4, title, ha="right", bold=True)
       a.plot([-2.2], [y], alpha=0)
       for i, (name, c) in enumerate(items):
           a.add_patch(Rectangle((i, y), 1, 0.8, fc=col[c], alpha=0.35, ec=C["grey"]))
           text(a, i + 0.5, y + 0.4, name)
           text(a, i + 0.5, y + 1.05, str(i), color=C["grey"])
   fig, (a1, a2) = strip(2, ["Vector of particles, and the cell map", "partition(cell ≤ c2)"], height=2.4)
   vec(a1, before, 1.6, "vector")
   for c, name in cells.items():
       idx = [str(i) for i, (_, cc) in enumerate(before) if cc == c]
       a1.add_patch(Rectangle((c * 1.5, 0), 1.4, 0.8, fc=col[c], alpha=0.2, ec=C["grey"]))
       text(a1, c * 1.5 + 0.7, 0.4, name + ": " + ",".join(idx))
   text(a1, -0.3, 0.4, "cell map", ha="right", bold=True)
   a1.plot([-2.2], [0], alpha=0)
   vec(a2, after, 0.5, "vector")
   brace(a2, 0, 4, 0.3, "returned range: pass", color=C["ok"], up=False)
   brace(a2, 4, 6, 0.3, "fail (c3)", color=C["bad"], up=False)
   finish(fig, [a1, a2])

* ``nbr_particles_in(box)``, ``export_particles(box, dest)`` only visit the
  cells of the box (:src:`src/core/data/particles/particle_array.hpp#L149-L179`).
* ``partition`` (:src:`src/core/utilities/cellmap.hpp#L400-L445`) is what the
  pusher's *selectors* use ("still in my domain?", "in the ghost layer?").
* ``push_back``, ``change_icell`` and ``erase`` on an index range keep the map
  up to date (:src:`src/core/data/particles/particle_array.hpp#L99-L190`); the
  iterator ``erase`` does not (comment at L101-L112).

The tricks
----------

Selections test cells, so the map must cover every reachable cell
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** A particle in a cell outside the map's box is never visited: a
selector would keep it by mistake. A half push can change a particle's cell by
up to 2.
**Fix.** Every array maps the ghost box grown by
``particleMapMargin`` = max(2, 4 − particle ghost width)
(:src:`src/amr/data/particles/particles_data.hpp#L68-L87`,
:src:`src/core/data/grid/gridlayout.hpp#L59-L68`). :pr:`1592`.

One SAMRAI variable for five arrays
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** Five variables per population would multiply the schedules.
**Fix.** One ``ParticlesData`` holds all five; the overlap type says what to
copy: ``ParticlesDomainOverlap`` = grab patch ghosts into my domain, cell
overlap = copy domain particles
(:src:`src/amr/data/particles/particles_data.hpp#L245-L270`); refinement goes
through the split operators (:doc:`/dev/particle_splitting`).

Restarts store four arrays out of five
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Domain + the three level ghost arrays are written; patch ghosts are empty
between steps (``putToRestart``,
:src:`src/amr/data/particles/particles_data.hpp#L114-L141`).

Go deeper
---------

* :doc:`/dev/grabbing_particles` · :doc:`/dev/level_ghost_particles` · :doc:`/dev/moments`

Related
-------

* :doc:`/dev/resources_manager` · :doc:`/dev/ppc_order` · :doc:`/numerics/particles`
