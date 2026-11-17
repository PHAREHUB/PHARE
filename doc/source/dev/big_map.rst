.. _dev-big-map:

===========
The big map
===========

.. container:: breadcrumb

   You are here: :doc:`Developer guide </dev/index>` › The big map

**In one sentence:** PHARE is a Python front end (input, run control, data
analysis) on top of a C++ engine in layers: plain physics and numerics at the
bottom, mesh refinement through the :term:`SAMRAI` library above it, and the
simulator that ties them together.

Each box below opens the page that explains that part.


The picture
-----------

Running a simulation, from the script down to the physics (click a box):

.. flow::

   Input script | pyphare.pharein | /dev/python_to_cpp
   fills
   PHAREDict | src/initializer, filled through the pybind module dictator | /dev/python_to_cpp
   then
   Python Simulator | pyphare.simulator: setup, advance | /dev/one_time_step
   calls, through pybind modules
   C++ Simulator<options> | src/python3, src/simulator | /dev/one_time_step
   uses
   AMR layer | src/amr, on top of SAMRAI: messengers, solvers, integrator, tagging | /dev/messengers
   uses
   Core layer | src/core: fields, particles, numerical schemes | /dev/gridlayout

Writing and reading the output:

.. flow::

   C++ Simulator | src/simulator
   writes
   Diagnostics, restarts | src/diagnostic, src/restarts, src/hdf5 | /dev/diagnostics_restarts
   HDF5 files |
   read by
   pharesee | pyphare.pharesee: read and plot | /data_analysis/reading_data

Nothing in a lower box calls a box above it: ``src/core`` is used by ``src/amr`` and uses nothing else of PHARE except
the dictionary.


How it works
------------

#. **You write a Python script** with ``pyphare.pharein``: a ``Simulation``,
   a model, diagnostics. Running it fills a C++ dictionary of parameters,
   the :term:`PHAREDict`. See :doc:`/dev/python_to_cpp`.
#. **The Python** ``Simulator`` **loads one compiled** (:term:`pybind`) **module** that matches the
   options of the run (dimension, interpolation order, ...), builds the
   SAMRAI hierarchy and the C++ simulator, then calls ``advance()`` in a loop.
   See :doc:`/dev/one_time_step`.
#. **The C++ simulator** (``src/simulator``) owns the mesh hierarchy, the
   time integrator, the diagnostics and the restarts.
#. **The amr layer** (``src/amr``) is where PHARE meets SAMRAI: it stores
   PHARE data on SAMRAI :term:`patches <patch>`, moves it between patches and
   :term:`levels <level>` (:term:`messengers <messenger>`), advances each level
   with a solver, and decides where to refine.
#. **The core layer** (``src/core``) holds the physics: field and particle
   containers, the grid layout, and the numerical schemes (pusher,
   deposition, Faraday, Ampère, Ohm, the MHD finite volume scheme). It works
   on one patch at a time and knows nothing about SAMRAI.
#. **Output** is written by ``src/diagnostic`` and ``src/restarts`` into
   HDF5 files, which ``pyphare.pharesee`` reads back for analysis.


Where the code is
-----------------

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Folder
     - What is in it
   * - ``pyphare/pyphare/pharein``
     - The Python input API (``Simulation``, models, diagnostics) and the code
       that fills the ``PHAREDict``.
   * - ``pyphare/pyphare/simulator``
     - The Python ``Simulator`` class that drives a run.
   * - ``pyphare/pyphare/pharesee``
     - Reading and plotting output files.
   * - ``src/python3``
     - The pybind modules: one per :term:`build permutation`
       (``cpp_<options>``), plus ``cpp_etc`` for everything not tied to the
       options.
   * - ``src/initializer``
     - The ``PHAREDict`` and the code that fills it from Python: the pybind module ``dictator``
       (:src:`src/initializer/dictator.cpp`), through which Python writes the dictionary.
   * - ``src/simulator``
     - ``Simulator<options>``, the C++ object behind ``advance()``.
   * - ``src/amr``
     - SAMRAI wrappers, patch data, messengers, solvers, tagging, load
       balancing, the multiphysics integrator.
   * - ``src/core``
     - Physics data structures and numerical schemes, independent of SAMRAI.
   * - ``src/diagnostic``, ``src/restarts``, ``src/hdf5``
     - Writing output and restart files.
   * - ``src/phare/phare.cpp``
     - A stand-alone executable that runs a Python input script without the
       Python ``Simulator``.


The tricks
----------

The physics layer does not know about SAMRAI
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** Mesh refinement code is large and hard to test. If every
numerical routine depended on it, each one would need a full SAMRAI hierarchy
to be tested or reused.

**Fix.** ``src/core`` never includes SAMRAI or ``src/amr`` headers: its
classes work on one grid layout and plain arrays. The ``src/amr`` layer
attaches these objects to SAMRAI patches and calls them patch by patch. So
the core schemes can be unit-tested on a single grid, and the same code runs
on every patch of every level.

The options of a run are fixed when PHARE is compiled
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** Dimension, interpolation order and the other main options
change the type of almost every object. Choosing them at run time with
``if`` statements inside hot loops would be slow.

**Fix.** They are template parameters of ``Simulator<options>``. Each
supported combination is compiled once into its own Python module, and the
Python side imports the one that matches the input script. See
:doc:`/dev/python_to_cpp`.


Go deeper
---------

* :doc:`/dev/one_time_step`: what happens in one call to ``advance()``.
* :doc:`/dev/python_to_cpp`: how the input script becomes a C++ simulation.
* :doc:`/dev/gridlayout`: where each quantity sits on a patch.
* :doc:`/dev/resources_manager`: how PHARE data lives on SAMRAI patches.
* :doc:`/dev/messengers`: how data moves between patches and levels.

Related
-------

* :doc:`/dev/working_on_phare`: building, testing, conventions.
* :doc:`/dev/cpp_api`: the generated C++ reference.
