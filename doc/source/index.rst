:orphan:

.. PHARE documentation master file, created by
   sphinx-quickstart on Sat Oct 30 17:25:13 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


PHARE
-----



PHARE is a Hybrid Particle-In-Cell (PIC) code.
It solves the evolution of the Vlasov equation of an arbitrary number of
ion populations in a Lagrangian way. Electrons are modeled as a single fluid.
Their momentum equation is used to compute the electric field, assuming quasineutrality.

Using Adaptive Mesh Refinement, provided by the library SAMRAI, PHARE aims at
filling the gap between sub-ion scales and large "MHD" scales by increasing
the mesh resolution wherever the solution needs it.



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


Theory
------
.. toctree::
   :caption: THEORY
   :maxdepth: 1
   :hidden:

   theory/hybridpic
   theory/pic
   theory/spatial_discretization
   theory/temporal_discretization
   theory/amr




Getting and Building
--------------------
.. toctree::
   :caption: BUILDING
   :maxdepth: 1
   :hidden:

   getting
   build



Usage
-----
.. toctree::
   :caption: USAGE
   :maxdepth: 1
   :hidden:

   usage/simulation_inputs
   usage/run_from_python
   usage/run_from_exe
   usage/examples




Data analysis
-------------
.. toctree::
   :caption: ANALYSIS
   :maxdepth: 1
   :hidden:

   pharesee/get_data
   pharesee/plotting_fields
   pharesee/plotting_distributions




Development
-----------
.. toctree::
   :caption: DEVELOPMENT
   :maxdepth: 1
   :hidden:

   development/tests

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


