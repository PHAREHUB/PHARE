:orphan:

.. PHARE documentation master file, created by
   sphinx-quickstart on Sat Oct 30 17:25:13 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


PHARE
-----



PHARE is a Hybrid Particle-In-Cell (PIC) code with MHD capabilities.
It solves the evolution of the Vlasov equation of an arbitrary number of
ion populations in a Lagrangian way. Electrons are modeled as a single fluid.
Their momentum equation is used to compute the electric field, assuming quasineutrality.

Using Adaptive Mesh Refinement, provided by the library SAMRAI, PHARE aims at
filling the gap between sub-ion scales and large "MHD" scales by increasing
the mesh resolution wherever the solution needs it.

PHARE also includes a finite volume MHD solver that can run standalone or
coupled with the Hybrid PIC solver on coarser AMR levels.


.. raw:: html

   <style>
   /* front page: hide chapter titles
    * needed for consistent HTML-PDF-EPUB chapters*/
   div#getting-and-building.section,
   div#usage.section,
   div#theory.section,
   div#data-analysis.section,
   div#development.section,
   div#indices-and-tables
   {
       display:none;
   }
   </style>




.. Theory
.. ------
.. toctree::
   :caption: THEORY
   :maxdepth: 1
   :hidden:

   theory/hybridpic
   theory/pic
   theory/spatial_discretization
   theory/temporal_discretization
   theory/amr
   theory/mhd




.. Getting and Building
.. --------------------
.. toctree::
   :caption: BUILDING
   :maxdepth: 1
   :hidden:

   getting
   build



.. Usage
.. -----
.. toctree::
   :caption: USAGE
   :maxdepth: 1
   :hidden:

   usage/simulation_inputs
   usage/run_from_python
   usage/run_from_exe



.. Tutorials
.. ---------
.. toctree::
   :caption: TUTORIALS
   :maxdepth: 1
   :hidden:

   tutorials/alfven_wave_1d
   tutorials/harris_2d
   tutorials/multi_population_1d
   tutorials/mhd_shock_1d
   tutorials/mhd_orszag_tang_2d
   tutorials/amr_tagging



.. Data analysis
.. -------------
.. toctree::
   :caption: ANALYSIS
   :maxdepth: 1
   :hidden:

   pharesee/get_data
   pharesee/plotting_fields
   pharesee/plotting_distributions
   pharesee/reconnection_analysis




.. Indices and tables
.. ------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

