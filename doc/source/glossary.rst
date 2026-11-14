========
Glossary
========

Short definitions of the words used in this documentation. Each entry says
what the word means *in PHARE*. Pages link here with :term:`patch` style
links; follow the "see" links for the full story.

Physics and models
==================

.. glossary::

   hybrid PIC
   hybrid model
      A plasma model where ions are followed as particles (particle-in-cell)
      and electrons are a massless fluid. The electric field comes from the
      electron momentum equation (:term:`Ohm's law`). See
      :doc:`/theory/hybridpic`.

   MHD
   magnetohydrodynamics
      A model where the whole plasma is one conducting fluid. PHARE solves it
      with a :term:`finite volume` method. See :doc:`/theory/mhd`.

   quasineutrality
      The assumption that electron and ion charge densities are equal at every
      point. The hybrid model uses it to get the electron density from the ion
      density.

   Ohm's law
      In the hybrid model, the electron momentum equation solved for the
      electric field. It has an ideal term with the electron velocity (which
      contains the Hall term), an electron pressure term, and optional
      resistive and hyper-resistive terms. See :doc:`/dev/electromagnetic`.

   closure
      The rule that gives the electron pressure, which the fluid equations
      cannot provide by themselves (for example isothermal electrons). See
      :doc:`/user_guide/electrons`.

   population
      A group of ions with the same charge and mass and its own initial
      distribution (for example "protons" or "beam"). Each population has its
      own particles and its own :term:`moments`.

   moments
      Density and flux (density times velocity) computed from the particles
      of a population, on the grid. The momentum tensor is computed too, for
      diagnostics only. See :doc:`/dev/moments`.

   normalization
      PHARE works in dimensionless units built from a reference magnetic field
      and density. See :doc:`/theory/normalization`.

Particles
=========

.. glossary::

   macro-particle
      One computational particle. It stands for many real ions, through its
      weight.

   interpolation order
      The order of the shape function that spreads a particle over the grid
      nodes (1, 2 or 3). Higher order means smoother, more expensive.
      See :doc:`/numerics/particles`.

   deposition
      Adding each particle's contribution to the density and flux on the grid
      nodes around it, with the particle's shape function.

   pusher
      The routine that moves particles one time step in the electric and
      magnetic fields (Boris scheme). See :doc:`/numerics/particles`.

   particle splitting
      Replacing one coarse particle by several smaller ones on a finer level,
      so that particle count per cell stays reasonable after
      refinement. See :doc:`/dev/particle_splitting`.

   refined particle number
      How many fine particles one coarse particle is split into
      (``refined_particle_nbr``). The values available depend on dimension and
      interpolation order.

   patch ghost particle
      A particle of a patch that has just moved out of the patch into cells
      of a same-level neighbour (the patch's ghost layer). The patch deposits
      it, then the neighbour takes it into its own domain, and the patch ghost
      array is emptied. See :doc:`/dev/grabbing_particles`.

   level ghost particle
      A particle in the ghost layer at the border of a refined level, made by
      splitting particles of the coarser level. See
      :doc:`/dev/level_ghost_particles`.

Mesh and adaptive mesh refinement
=================================

.. glossary::

   AMR
   adaptive mesh refinement
      Using a finer mesh only where the solution needs it. The fine mesh is
      made of :term:`patches <patch>` grouped in :term:`levels <level>`.
      See :doc:`/numerics/amr`.

   SAMRAI
      The C++ library (from Lawrence Livermore National Laboratory) that
      PHARE uses to manage patches, levels, communication between them and
      the recursive time stepping.

   level
      All the patches with the same mesh size. Level 0 covers the whole
      domain; each finer level is 2 times finer (the :term:`refinement
      ratio`) and covers only part of the level below.

   patch
      A rectangular block of cells on one level. Each MPI process owns some
      patches. A patch is the unit of data and of work.

   patch hierarchy
      The set of all levels and their patches at a given time. Output files
      are organized the same way. See :doc:`/data_analysis/patch_hierarchy`.

   refinement ratio
      Ratio of mesh sizes between two consecutive levels. It is 2 in PHARE
      (:src:`src/amr/amr_constants.hpp`).

   ghost cell
      A cell just outside a patch that holds a copy or an interpolation of
      data owned by a neighbour patch or by the coarser level, so that the patch can compute derivatives and move particles near its
      edge. See :doc:`/dev/ghost_cells`.

   patch ghost
      A ghost cell filled from a neighbour patch on the same level.

   level ghost
      A ghost cell at the border of a refined level, filled from the coarser
      level by refinement (interpolation in space and time).

   refinement
      Building fine-level values from the coarser level (fields:
      interpolation; particles: :term:`particle splitting`). See
      :doc:`/dev/field_refinement`.

   coarsening
      Replacing coarse-level values by values taken from the finer level, where
      a finer level exists. In the hybrid solver the electric field is averaged
      along each edge from the fine edges, the ion moments are copied from the
      fine node that sits on the coarse node, and the magnetic field is not
      coarsened: it is recomputed by :term:`reflux`. See :doc:`/dev/coarsening`.

   tagging
      Marking the cells that need a finer mesh, from a criterion on the
      solution. See :doc:`/dev/tagging_regridding`.

   sub-cycling
      A finer level advancing with several smaller time steps for one step of
      the level below. PHARE asks SAMRAI for fine steps no larger than the
      coarse step divided by 4 (refinement ratio 2, squared); with a fixed
      coarse step, always the case in hybrid, this gives exactly 4 sub-steps.
      See :doc:`/dev/subcycling`.

   time interpolation
      Setting a fine-level ghost value at a time between two coarse-level
      times, as a linear mix :math:`(1-\alpha)\,\text{old} + \alpha\,\text{new}`:
      *old* is saved by the coarse level at the start of its step, *new* is
      the coarse value at the end of it. See :doc:`/dev/time_interpolation`.

   synchronization
      The moment a finer level has caught up in time with the level below,
      which then takes the finer level's result. Hybrid coarsens E and the ion
      moments and recomputes B by :term:`reflux`; MHD redoes its coarse step
      with the fine fluxes. See :doc:`/dev/refluxing` and :doc:`/dev/coarsening`.

   regridding
      Rebuilding the finer levels from the tagged cells, so that refined
      regions follow the solution.

   fill pattern
      A rule that decides which values a data copy between patches may write.
      See :doc:`/dev/fill_patterns`.

   schedule
      A SAMRAI object that knows which pieces of data go from which patch to
      which patch, prepared once and run many times.

   messenger
      The PHARE object that owns the schedules and runs them at the right
      moment of a time step (filling ghosts, refining, coarsening). See
      :doc:`/dev/messengers`.

   reflux
      A correction made when a fine level has caught up in time with the
      coarser level: the coarse update is redone with the fine-level fluxes,
      summed over the fine sub-steps, so that both levels agree. The MHD
      solver uses the fluid fluxes; the hybrid solver uses the electric field,
      which plays the role of the flux of :math:`B` in Faraday's law.

Grid and numerical methods
==========================

.. glossary::

   Yee layout
   staggered grid
      Placing different field components at different points of a cell
      (nodes, faces, edges, centers), so that discrete curls are
      second-order accurate and the discrete divergence of a curl is exactly
      zero. See :doc:`/numerics/grid`.

   primal
   dual
      Along one direction, a quantity is *primal* if it sits on cell
      boundaries (nodes) and
      *dual* if it sits on cell centers. See :doc:`/dev/gridlayout`.

   GridLayout
      The C++ class that knows where each quantity sits on the grid and how to
      convert between index spaces. See :doc:`/dev/gridlayout`.

   predictor-predictor-corrector
   PPC
      The time scheme of the hybrid solver: two prediction stages and one
      correction stage per time step. See :doc:`/numerics/time_integration`
      and :doc:`/dev/ppc_order`.

   finite volume
      A method that evolves cell averages using fluxes through the cell
      faces. Used by the MHD solver. See :doc:`/numerics/mhd_finite_volume`.

   Godunov method
      A finite volume method where the flux through each face comes from
      solving a :term:`Riemann problem` between the two neighbouring cells.

   Riemann problem
   Riemann solver
      The evolution of two uniform states in contact; a Riemann solver gives
      (often approximately) the flux it produces at the interface.

   constrained transport
   CT
      Updating the magnetic field from electric fields on cell edges. The
      update does not change the discrete divergence of :math:`B`: if it is
      zero at the start, it stays zero to round-off.

   CFL condition
      The time step limit that keeps an explicit scheme stable: information
      must not cross more than about one cell per time step.

Running and code
================

.. glossary::

   MPI
      The standard library for running one program as many processes
      (*ranks*) that exchange messages. PHARE runs in parallel with MPI only.

   rank
      One MPI process.

   collective call
      An MPI operation that every rank must make, in the same order. If one
      rank skips it, the others wait forever. PHARE's helpers for these are
      in :src:`src/mpi/mpi_utils.hpp#L22-L47` (for example ``mpi::collect``).

   load balancing
      Moving patches between ranks so that each rank has a similar amount of
      work. See :doc:`/user_guide/load_balancing`.

   diagnostics
      Output files (HDF5) written during a run. See
      :doc:`/user_guide/diagnostics`.

   restart
      Starting a run again from files written by an earlier run. See
      :doc:`/user_guide/restarts`.

   permutation
   build permutation
      One compiled combination of compile-time options (dimension,
      interpolation order, refined particle number, MHD scheme). Each
      permutation is a separate Python module, listed in
      :src:`res/sim/all.txt`. See :doc:`/dev/python_to_cpp`.

   pybind
      pybind11, the library that makes C++ classes callable from Python.

   PHAREDict
      The C++ dictionary that holds every input parameter, filled from the
      Python input script. See :doc:`/dev/python_to_cpp`.
