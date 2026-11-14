===============
Developer guide
===============

How PHARE works inside, for people who read and change its C++ and Python
sources. Every page of this guide is linked from this one.

To *run* simulations, start from the :doc:`user documentation </index>`
instead.

How to use this guide
=====================

- **Team member refreshing one trick:** find its group below and go straight
  to the page.
- **Newcomer:** read the pages of the *Start* box in order, then follow the
  links at the bottom of each page.
- Each page links down to finer pages (*Go deeper*) and back up to the page
  above it (the line at the top).
- Code links point to the exact lines at the commit the pages were checked
  against.

.. container:: dev-dashboard

   .. container:: dev-card

      **Start**

      #. :doc:`big_map`
      #. :doc:`one_time_step`
      #. :doc:`python_to_cpp`

   .. container:: dev-card

      **Data structures**

      - :doc:`gridlayout`
      - :doc:`fields`
      - :doc:`resources_manager`
      - :doc:`particles_ions`

   .. container:: dev-card

      **Ghosts and mesh refinement**

      - :doc:`ghost_cells`
      - :doc:`fill_patterns`
      - :doc:`field_refinement`
      - :doc:`coarsening`
      - :doc:`messengers`
      - :doc:`tagging_regridding`

   .. container:: dev-card

      **Particles**

      - :doc:`level_ghost_particles`
      - :doc:`particle_splitting`
      - :doc:`grabbing_particles`
      - :doc:`moments`
      - :doc:`ppc_order`

   .. container:: dev-card

      **Solvers**

      - :doc:`electromagnetic`
      - :doc:`mhd_solver`

   .. container:: dev-card

      **Input/output and performance**

      - :doc:`diagnostics_restarts`
      - :doc:`load_balancing`
      - :doc:`hpc_internals`

   .. container:: dev-card

      **Working on PHARE**

      - :doc:`working_on_phare`
      - :doc:`how_to_add`
      - :doc:`cpp_api`

.. The toctree below is hidden: it only tells Sphinx that these pages belong
   to the developer guide (previous/next links, no "not in any toctree"
   warning). The left menu shows a single "Developer guide" entry.

.. toctree::
   :hidden:

   big_map
   one_time_step
   python_to_cpp
   gridlayout
   fields
   resources_manager
   particles_ions
   ghost_cells
   fill_patterns
   field_refinement
   coarsening
   messengers
   tagging_regridding
   level_ghost_particles
   particle_splitting
   grabbing_particles
   moments
   ppc_order
   electromagnetic
   mhd_solver
   diagnostics_restarts
   load_balancing
   hpc_internals
   working_on_phare
   how_to_add
   cpp_api
