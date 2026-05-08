
==============
Building PHARE
==============


Production Build
----------------

.. code-block:: bash

    cd path/to/dir/containing/PHARE
    mkdir build && cd build
    cmake -DCMAKE_BUILD_TYPE=Release \
          -DCMAKE_CXX_FLAGS="-O3 -march=native -mtune=native" \
          ../PHARE
    make -j$(nproc)


Debug Build
-----------

.. code-block:: bash

    cd path/to/dir/containing/PHARE
    mkdir build && cd build
    cmake -DCMAKE_BUILD_TYPE=Debug \
          -DCMAKE_CXX_FLAGS="-g3 -O0 -DPHARE_DIAG_DOUBLES=1" \
          ../PHARE
    make -j$(nproc)


CMake Options
-------------

.. list-table::
   :header-rows: 1
   :widths: 28 12 60

   * - Name
     - Default
     - Description
   * - ``devMode``
     - OFF
     - Enable strict compiler flags (``-Werror``); recommended during development
   * - ``test``
     - ON
     - Build the GoogleTest C++ test suite
   * - ``testMPI``
     - OFF
     - Run tests through ``mpirun`` instead of directly
   * - ``PHARE_MPI_PROCS``
     - 2
     - Number of MPI processes used when ``testMPI=ON``
   * - ``lowResourceTests``
     - OFF
     - Skip computationally heavy tests; useful in CI or on laptops
   * - ``asan``
     - OFF
     - Instrument with AddressSanitizer
   * - ``ubsan``
     - OFF
     - Instrument with UndefinedBehaviorSanitizer
   * - ``coverage``
     - OFF
     - Generate coverage reports via gcovr
   * - ``bench``
     - OFF
     - Build micro-benchmarks
   * - ``SAMRAI_ROOT``
     - (auto)
     - Path to a pre-installed SAMRAI; skips the automatic download
   * - ``PHARE_EXEC_LEVEL_MIN``
     - 1
     - Minimum test execution level (1–10)
   * - ``PHARE_EXEC_LEVEL_MAX``
     - 10
     - Maximum test execution level (1–10)

Pass any option with ``-D``, for example:

.. code-block:: bash

    cmake -Dtest=ON -DdevMode=ON -DtestMPI=ON -DPHARE_MPI_PROCS=4 ../PHARE


Setting up the Environment
--------------------------

After a successful build, the Python bindings live inside the build tree. You must add
both the ``pyphare`` package directory and the build directory to ``PYTHONPATH`` so that
Python can find them.

.. code-block:: bash

    export PYTHONPATH=/path/to/PHARE/pyphare:/path/to/build:$PYTHONPATH

.. warning::

   Forgetting to set ``PYTHONPATH`` is the most common reason ``import pyphare`` fails.
   Add the export to your shell profile or to a module file on HPC systems so it is
   always present when running simulations.


Verifying the Build
-------------------

.. code-block:: bash

    python3 -c "import pyphare.pharein as ph; print('PHARE ready')"

A successful output is simply::

    PHARE ready

You can also run the full test suite with CTest:

.. code-block:: bash

    cd /path/to/build
    ctest -j4 --output-on-failure

To run a single test by name:

.. code-block:: bash

    ctest -R test_name -V


Troubleshooting
---------------

HDF5 is Not the Parallel Build
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Symptom**: Link errors mentioning missing ``H5_HAVE_PARALLEL`` or failures at runtime
when writing diagnostic files.

**Fix**: Install the MPI-aware HDF5 variant and point CMake to it:

.. code-block:: bash

    # Debian/Ubuntu example
    sudo apt install libhdf5-openmpi-dev
    cmake -DHDF5_ROOT=/usr/lib/x86_64-linux-gnu/hdf5/openmpi ../PHARE

Wrong Python Interpreter
~~~~~~~~~~~~~~~~~~~~~~~~

**Symptom**: ``import pyphare`` succeeds but the C++ extension module is missing, or
Python picks up a system interpreter that differs from the one used at configure time.

**Fix**: Activate your virtual environment before configuring **and** before running, and
pin the interpreter explicitly:

.. code-block:: bash

    cmake -DPython3_EXECUTABLE=$(which python3) ../PHARE

MPI Not Found
~~~~~~~~~~~~~

**Symptom**: CMake reports ``Could NOT find MPI``.

**Fix**: Ensure ``mpicc`` / ``mpicxx`` are on ``PATH``. On systems with environment
modules:

.. code-block:: bash

    module load mpi/openmpi-x86_64
    cmake ../PHARE

Alternatively, set ``MPI_HOME``:

.. code-block:: bash

    cmake -DMPI_HOME=/path/to/mpi ../PHARE

SAMRAI Download Fails
~~~~~~~~~~~~~~~~~~~~~

**Symptom**: Configure step stops with a network error while fetching SAMRAI.

**Fix**: Download and build SAMRAI manually on a machine with internet access, then pass
its installation prefix at configure time:

.. code-block:: bash

    cmake -DSAMRAI_ROOT=/path/to/samrai-install ../PHARE

.. note::

   On HPC clusters it is common to pre-build SAMRAI once and reuse it across multiple
   PHARE builds. This also speeds up configure time significantly.
