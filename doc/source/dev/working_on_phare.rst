.. _dev-working_on_phare:

================
Working on PHARE
================

.. container:: breadcrumb

   You are here: :doc:`/dev/index` › Working on PHARE.
   Installing for users: :doc:`/getting_started/install`.

**In one sentence:** how to configure a developer build, run the tests you
need, and get a pull request through the checks.

The picture
-----------

.. mermaid::

   flowchart TD
     I["issue"] --> B["branch + code + tests"]
     B --> C["cmake configure<br/>(devMode, test levels)"]
     C --> M["build"] --> T["ctest<br/>(-R / -L / -LE)"]
     T --> PR["pull request<br/>(template)"]
     PR --> G["gate: diff checks"]
     PR --> CI["builds + tests<br/>(GitHub Actions, TeamCity)"]
     G --> R["review"]
     CI --> R

Configure a developer build
---------------------------

PHARE is built with CMake. MPI and parallel HDF5 must be installed
(``-DHDF5_ROOT=...`` points to HDF5). Other dependencies (SAMRAI, pybind11,
HighFive, cppdict, googletest) are cloned into ``subprojects/`` at configure
time unless you point CMake to an installed copy (``-DSAMRAI_ROOT=...``,
``-DGTEST_ROOT=...``).

.. code-block:: bash

   mkdir -p build/dev && cd build/dev
   cmake ../.. -DCMAKE_BUILD_TYPE=Debug -DdevMode=ON \
               -DCMAKE_CXX_FLAGS="-O3 -DPHARE_DIAG_DOUBLES=1"
   cmake --build . -j 8

This is close to what the GitHub CI build uses
(:src:`.github/workflows/cmake_ubuntu.yml#L64-L71`).

.. list-table:: Options developers use most (:src:`res/cmake/options.cmake`)
   :header-rows: 1
   :widths: 35 15 50

   * - Option
     - Default
     - Effect
   * - ``CMAKE_BUILD_TYPE``
     - ``Release``
     - ``Debug`` keeps PHARE's debug checks (``PHARE_DEBUG_DO``, compiled
       when ``NDEBUG`` is not defined) and SAMRAI's assertions when SAMRAI
       is built in-tree. PHARE adds no ``-O`` or ``-g`` itself; CMake's
       defaults apply (with gcc/clang, Release: ``-O3 -DNDEBUG``, Debug:
       ``-g``). CI builds Debug with ``-O3`` to keep the checks and still
       run fast.
   * - ``devMode``
     - OFF
     - Warnings become errors (:src:`res/cmake/def.cmake#L74-L80`), errors
       become exceptions instead of ``MPI_Abort`` (see :doc:`hpc_internals`,
       :src:`res/cmake/def.cmake#L25-L27`), ccache if found. Also, every
       configure runs ``git pull`` in each cloned subproject
       (:src:`res/cmake/dep.cmake#L10-L13`): SAMRAI can change under you.
   * - ``test`` / ``testMPI``
     - ON / OFF
     - Register tests. With ``testMPI``, C++ tests run under ``mpirun``
       with ``PHARE_MPI_PROCS`` ranks (default 2,
       :src:`res/cmake/def.cmake#L179-L184`), serial-only tests are skipped,
       and Python runs registered with ``phare_python3_exec`` use 2 ranks.
   * - ``PHARE_EXEC_LEVEL_MIN`` / ``_MAX``
     - 1 / 10
     - Which tests are registered (next section).
   * - ``asan`` / ``ubsan``
     - OFF
     - Address / undefined-behaviour sanitizers.
   * - ``withPhlop`` / ``withCaliper``
     - OFF
     - Timer backends, see :doc:`/hpc/measuring_performance`.
   * - ``bench``
     - OFF
     - Build the micro-benchmarks (google-benchmark).
   * - ``-DPHARE_DIAG_DOUBLES=1`` (in ``CMAKE_CXX_FLAGS``)
     - single precision
     - Write diagnostics (both formats) in double instead of single
       precision (:src:`src/diagnostic/detail/h5writer.hpp#L38`,
       :src:`src/diagnostic/detail/vtkh5_type_writer.hpp#L137`).

Tests
-----

Tests are registered with CTest. List them with ``ctest -N``.

Kinds of tests
   ``tests/core`` and ``tests/amr``: C++ unit tests (googletest).
   ``tests/simulator``: Python-driven runs of the whole simulator (restarts,
   diagnostics, load balancing...). ``tests/functional``: physics runs
   (dispersion, Harris, shock, MHD cases...). ``pyphare/pyphare_tests``:
   Python unit tests.

Levels: which tests exist in your build
   Heavy tests are registered with a *level*; a test exists only if
   ``PHARE_EXEC_LEVEL_MIN <= level <= PHARE_EXEC_LEVEL_MAX``
   (``phare_exec`` / ``phare_mpi_python3_exec``,
   :src:`res/cmake/def.cmake#L318-L343`). Plain unit tests count as level
   1. The defaults (1 to 10) give unit tests, simulator tests (level 9) and
   the hybrid dispersion tests (level 10). Heavier functional runs use 11,
   21, 101 or 121: configure with ``-DPHARE_EXEC_LEVEL_MAX=100`` to get
   levels up to 100, ``-DPHARE_EXEC_LEVEL_MAX=121`` (or more) to get them
   all.

Labels: which tests you run
   A label tags tests inside a build. At the pinned commit there are two:
   ``physics`` (dispersion tests) and ``nightly`` (3D particle
   initialization), set with ``phare_test_labels``
   (:src:`res/cmake/def.cmake#L352-L360`).

.. code-block:: bash

   ctest -j 8 --output-on-failure        # everything registered
   ctest -R restarts                     # tests whose name matches
   ctest -L physics                      # only the physics label
   ctest -LE "nightly|physics"           # all but those labels

Check scripts without running them
   ``PHARE_DRY_RUN=1`` builds and validates each simulation but makes
   ``initialize``, ``advance`` and ``run`` do nothing
   (:src:`pyphare/pyphare/pharein/simulation.py#L795-L797`), so ``ctest``
   only checks that scripts set up correctly. CI uses it on the heavy levels
   (:src:`res/ci/jobs/pull_request/dry_run_functional_tests.sh#L10-L12`).

Adding a test
   In the test folder's ``CMakeLists.txt``, use ``add_phare_test`` (C++),
   ``phare_python3_exec(level name script.py dir)`` or
   ``phare_mpi_python3_exec(level N name script.py dir)`` for a Python run
   on ``N`` ranks. For such a Python run, ``phare_test_labels(name label)``
   makes it selectable by label (it only matches Python tests).

.. admonition:: In preparation
   :class: in-preparation

   Typical run time of each test level on a workstation.

Continuous integration
----------------------

Every pull request goes through:

* **Gate** (``.github/workflows/gate.yml``): checks on the diff only, no
  build: skipped tests, loosened test tolerances, changed input keywords,
  formatting (clang-format), Python lint (ruff), removed tests, the PR
  description, cyclic imports, clang-tidy, copy-paste, workflow files,
  dependencies. Its scripts are
  taken from ``master``, so a PR cannot change the rules it is judged by
  (:src:`.github/workflows/gate.yml#L1-L14`).
* **Builds and tests**: a GitHub Actions build with gcc and clang running
  ``ctest`` (``cmake_ubuntu.yml``), a coverage build (``coverage.yml``,
  ``ctest -LE "nightly|physics"``), and builds on the project's TeamCity
  server (configuration not stored in the repository).
* **Review threads**: a check fails while a review thread started by a
  person (not a bot) is unresolved (``review-threads.yml``).

Conventions
-----------

* Contribution rules: an issue first, one issue per PR, small PRs, tests
  with the code — see ``CONTRIBUTING.md`` and the PR template
  (``.github/pull_request_template.md``: Issue / What this implements /
  Design / Tests).
* C++ formatting: ``.clang-format`` at the root (clang-format version
  pinned in ``pyproject.toml``, ``lint`` group).
* Python: ``ruff.toml``; Python dependencies in ``pyproject.toml``; do not
  import ``h5py``, ``mpi4py`` or ``scipy.optimize`` at the top of a module
  (they load system libraries that may not match the build). See
  ``doc/conventions.md``.

Go deeper
---------

* :doc:`how_to_add` — recipes for common changes.
* :doc:`cpp_api`

Related
-------

* :doc:`/getting_started/build` — user build instructions.
* :doc:`/hpc/measuring_performance`
