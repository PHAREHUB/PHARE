===========
Build PHARE
===========

**Goal:** compile PHARE's C++ core and check that Python finds it.

You need the tools and the source from :doc:`install`.


Minimal build
-------------

From the folder that contains ``PHARE``:

.. code-block:: bash

   cmake -S PHARE -B build -DCMAKE_BUILD_TYPE=Release -Dtest=OFF
   cmake --build build -j 8

* ``-S PHARE -B build``: the source is in ``PHARE``, the build goes in a new
  folder ``build``.
* ``-DCMAKE_BUILD_TYPE=Release``: optimized code, for real runs.
* ``-Dtest=OFF``: skip PHARE's own test suite (it is on by default and makes
  the build much longer).
* ``-j 8``: compile 8 files at a time; use about the number of cores of your
  machine.

The first ``cmake`` call downloads SAMRAI, HighFive, cppdict and (if not
installed) pybind11 into
``PHARE/subprojects`` (internet needed), and SAMRAI is compiled together with
PHARE. Expect the build to take a long time: it compiles one version of the
core per kind of simulation (next section).

Then tell Python where PHARE is (both folders are needed):

.. code-block:: bash

   export PYTHONPATH=/path/to/PHARE/pyphare:/path/to/build:$PYTHONPATH

This line only holds in the current terminal. Put it in your shell start-up
file (for example ``~/.bashrc``) to keep it.

**Check:** this prints the kinds of simulation your build can run.

.. code-block:: bash

   python3 -c "from pyphare import cpp; print(cpp.built_simulator_ids())"

You should see a list like ``['1_1_2', '1_1_3', '1_2_2', ...]``. If it prints
``None``, Python does not find the build folder: check ``PYTHONPATH``. Now go to
:doc:`first_simulation`.


What the build makes
--------------------

For speed, PHARE compiles a separate version of its core for each combination
of a few settings that must be known at compile time: the number of dimensions,
and, for hybrid runs, the interpolation order and the number of particles a
particle is split into when it enters a finer mesh level; for :term:`MHD`
runs, the numerical scheme (reconstruction, slope limiter, Riemann solver,
Hall term).
Each combination is one Python module in ``build/pybindlibs``. When you run a
script, PHARE reads your ``Simulation`` settings and loads the matching module.

.. mermaid::

   flowchart TD
     F["res/sim/all.txt<br/>one line per combination<br/>e.g. 1,1,2"] --> B["build"]
     B --> M["build/pybindlibs/<br/>cpp_1_1_2, cpp_2_1_4, ...<br/>one module per line"]
     S["your script<br/>Simulation(cells=200, interp_order=1, ...)"] --> P["pick module<br/>cpp_1_1_2"]
     M --> P

The list of combinations is the file ``PHARE/res/sim/all.txt``. Its lines are
``dim,interp_order,refined_particle_nbr`` for hybrid runs and
``dim,reconstruction,limiter,riemann,hall`` for MHD runs. The default list
covers 1D, 2D and 3D, and every combination used by PHARE's tests.

To build faster, or to add a combination, give your own file:

.. code-block:: bash

   cmake -S PHARE -B build -DCMAKE_BUILD_TYPE=Release -Dtest=OFF \
         -DPHARE_PERMUTATIONS=/path/to/my_list.txt

For example, a file with the single line ``1,1,2`` builds only what
:doc:`first_simulation` needs: on a 12-core workstation this build takes
about 2 minutes with ``-j 8``, SAMRAI included.


Build options
-------------

.. list-table::
   :header-rows: 1
   :widths: 30 50 20

   * - Option
     - Meaning
     - Default
   * - ``CMAKE_BUILD_TYPE``
     - ``Release`` for runs; ``Debug`` to find bugs (slow, with checks)
     - ``Release``
   * - ``test``
     - also build PHARE's test suite
     - ``ON``
   * - ``PHARE_PERMUTATIONS``
     - file listing the combinations to compile (see above)
     - ``res/sim/all.txt``
   * - ``SAMRAI_ROOT``
     - use a SAMRAI already installed in this folder instead of downloading it
     - not set
   * - ``HDF5_ROOT``
     - where to find HDF5 if CMake does not find it alone
     - not set
   * - ``PHARE_EXE``
     - also build a stand-alone executable, ``phare-exe``
     - ``OFF``
   * - ``-DPHARE_DIAG_DOUBLES=1`` inside ``CMAKE_CXX_FLAGS``
     - write output files in double precision (about 16 significant digits)
       instead of single precision (about 7)
     - single precision

Options are given to the first ``cmake`` call as ``-D<option>=<value>``.

With your own compiler flags, for example optimized for your processor and
with double-precision output, put all flags in **one**
``CMAKE_CXX_FLAGS`` (if you give it twice, only the last one counts):

.. code-block:: bash

   cmake -S PHARE -B build -DCMAKE_BUILD_TYPE=Release -Dtest=OFF \
         -DCMAKE_CXX_FLAGS="-march=native -DPHARE_DIAG_DOUBLES=1"

A debug build, to run PHARE in a debugger:

.. code-block:: bash

   cmake -S PHARE -B build_debug -DCMAKE_BUILD_TYPE=Debug -Dtest=OFF \
         -DCMAKE_CXX_FLAGS="-g3 -O0 -DPHARE_DIAG_DOUBLES=1"
   cmake --build build_debug -j 8


Common mistakes
---------------

* ``ModuleNotFoundError: No module named 'pybindlibs.cpp_3_1_7'``, followed by
  ``This run needs the build permutation``, the line ``3,1,7`` and the list of
  what is built
  — your script asks for a combination your build does not have. Either change
  your settings to one in the list, or add the line (here ``3,1,7``) to your
  combinations file and rebuild.
* ``ModuleNotFoundError: No module named 'pybindlibs'`` — Python does not find
  the build folder. Add it to ``PYTHONPATH``.
* ``No module named 'pybindlibs.cpp_...'`` although the build has it — you run
  a different Python from the one used to build. Activate the same virtual
  environment, or rebuild with the Python you use.
* The first ``cmake`` call fails while cloning SAMRAI, HighFive or pybind11 —
  no internet access. Run that step on a machine with access, or install SAMRAI
  yourself and pass ``-DSAMRAI_ROOT``.
* ``HDF5 NOT PARALLEL!`` when running with ``mpirun`` — PHARE was built with
  a serial HDF5. Install or load a parallel HDF5 and rebuild.
* CMake does not find HDF5 or MPI — load them first (on a cluster:
  ``module load ...``), or give ``-DHDF5_ROOT=...``.
