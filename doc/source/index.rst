:orphan:

PHARE
=====

PHARE is a plasma simulation code with two models:

- **hybrid particle-in-cell (PIC)**: ions are particles, electrons are a fluid.
  It solves the Vlasov equation for any number of ion populations by following
  macro-particles. The electron momentum equation gives the electric field,
  assuming quasineutrality.
- **magnetohydrodynamics (MHD)**: the whole plasma is one fluid.

Using adaptive mesh refinement (AMR), provided by the SAMRAI library, PHARE
refines the mesh only where the solution needs it. A single run can then
resolve sub-ion scales in small regions of a much larger domain.

.. note::

   This documentation is a work in progress. It is not complete and may
   contain errors.

Where to start
--------------

- **New to PHARE?** :doc:`getting_started/what_is_phare`, then
  :doc:`getting_started/install` and :doc:`getting_started/first_simulation`.
- **Setting up a run?** The *User guide* in the menu lists every input block.
- **Working on the code?** The :doc:`dev/index` is a separate part of this site,
  for people who read and change the C++ and Python sources.


.. toctree::
   :caption: Getting started
   :maxdepth: 1
   :hidden:

   getting_started/what_is_phare
   getting_started/install
   getting_started/build
   getting_started/first_simulation
   getting_started/first_result


.. toctree::
   :caption: User guide
   :maxdepth: 1
   :hidden:

   user_guide/simulation
   user_guide/choosing_a_model
   user_guide/hybrid_inputs
   user_guide/electrons
   user_guide/mhd_inputs
   user_guide/boundary_conditions
   user_guide/amr
   user_guide/diagnostics
   user_guide/running
   user_guide/restarts
   user_guide/load_balancing
   user_guide/environment_variables


.. toctree::
   :caption: Data analysis
   :maxdepth: 1
   :hidden:

   data_analysis/reading_data
   data_analysis/patch_hierarchy
   data_analysis/plotting_fields
   data_analysis/distributions


.. toctree::
   :caption: Examples
   :maxdepth: 1
   :hidden:

   examples/index


.. toctree::
   :caption: Theory
   :maxdepth: 1
   :hidden:

   theory/hybridpic
   theory/pic
   theory/mhd
   theory/normalization


.. toctree::
   :caption: Numerics
   :maxdepth: 1
   :hidden:

   numerics/grid
   numerics/time_integration
   numerics/particles
   numerics/mhd_finite_volume
   numerics/domain_boundary_conditions
   numerics/amr


.. toctree::
   :caption: Performance & HPC
   :maxdepth: 1
   :hidden:

   hpc/running_at_scale
   hpc/memory
   hpc/measuring_performance
   hpc/scaling_results


.. toctree::
   :caption: Validation
   :maxdepth: 1
   :hidden:

   validation/hybrid_dispersion


.. toctree::
   :caption: Reference
   :maxdepth: 1
   :hidden:

   glossary


.. toctree::
   :caption: For developers
   :maxdepth: 1
   :hidden:

   dev/index


Indices
-------

* :ref:`genindex`
* :ref:`search`
