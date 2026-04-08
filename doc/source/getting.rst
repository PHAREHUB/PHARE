
=============
Getting PHARE
=============

Prerequisites
-------------

The following tools and libraries must be installed on your system before building PHARE.

.. list-table::
   :header-rows: 1
   :widths: 30 30 40

   * - Dependency
     - Minimum version
     - Notes
   * - C++ compiler
     - GCC >= 10, Clang >= 13
     - Must support C++20
   * - CMake
     - 3.20.1
     - Build system
   * - MPI
     - Any recent release
     - OpenMPI, MPICH, or Intel MPI
   * - HDF5
     - Any recent release
     - **Must** be the parallel (MPI-enabled) build
   * - Python
     - 3.8
     - Required for bindings and test suite

.. warning::

   HDF5 must be compiled with MPI support (parallel HDF5). The serial HDF5 library will
   cause link errors. On most distributions the parallel variant is provided by a separate
   package (e.g. ``libhdf5-openmpi-dev`` on Debian/Ubuntu).


Auto-downloaded Dependencies
----------------------------

The following libraries are fetched automatically by CMake at configure time and do not
need to be installed manually.

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Library
     - Purpose
   * - SAMRAI
     - Structured Adaptive Mesh Refinement Application Infrastructure
   * - pybind11
     - Python / C++ binding layer
   * - HighFive
     - Header-only C++ HDF5 wrapper
   * - GoogleTest
     - C++ unit-test framework

.. note::

   An internet connection is required during the first ``cmake`` invocation. If you are on
   a cluster without outbound access, pre-install SAMRAI and point CMake to it with
   ``-DSAMRAI_ROOT=/path/to/samrai``. See :doc:`build` for details.


Getting the Source
------------------

Clone the repository together with all its submodules in one step:

.. code-block:: bash

    git clone --recursive https://github.com/PHAREHUB/PHARE

If you already cloned without ``--recursive``, initialise submodules afterwards:

.. code-block:: bash

    cd PHARE
    git submodule update --init --recursive


Python Dependencies
-------------------

It is strongly recommended to use a dedicated virtual environment.

.. code-block:: bash

    python3 -m venv .venv
    source .venv/bin/activate
    pip install -r requirements.txt

Key packages installed by ``requirements.txt``:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Package
     - Purpose
   * - numpy
     - Array operations used throughout the analysis layer
   * - scipy
     - Signal processing and interpolation utilities
   * - h5py
     - Python interface to HDF5 diagnostic output
   * - matplotlib
     - Plotting utilities for post-processing scripts
   * - dill
     - Extended pickling used for Python-defined initial conditions
   * - pyyaml
     - YAML configuration parsing

Once the source is on disk, proceed to :doc:`build`.
