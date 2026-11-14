.. _dev-working_on_phare:

================
Working on PHARE
================

.. container:: breadcrumb

   You are here: :doc:`/dev/index` › Working on PHARE. Installing for users:
   :doc:`/getting_started/install`.

**In one sentence:** configure a developer build, run the tests you need, and get a pull request
through the checks.

The usual path: an issue → a branch with code and tests → configure, build, ``ctest`` → a pull
request (template) → gate checks, builds and tests → review.

Configure a developer build
---------------------------

.. code-block:: bash

   mkdir -p build/dev && cd build/dev
   cmake ../.. -DCMAKE_BUILD_TYPE=Debug -DdevMode=ON \
               -DCMAKE_CXX_FLAGS="-O3 -DPHARE_DIAG_DOUBLES=1"
   cmake --build . -j 8

Close to the GitHub CI build (:src:`.github/workflows/cmake_ubuntu.yml#L64-L71`).

* **Dependencies.** MPI and parallel HDF5 must be installed (``-DHDF5_ROOT``). SAMRAI and googletest are
  cloned into ``subprojects/`` unless ``-DSAMRAI_ROOT`` / ``-DGTEST_ROOT`` point to an install; pybind11 is
  cloned if none is found (or with ``-DforceGetPybind=ON``); HighFive and cppdict are always cloned
  (:src:`res/cmake/dep/pybind.cmake#L24-L34`, :src:`res/cmake/dep/cppdict.cmake#L5-L7`).

.. list-table:: Options developers use most (:src:`res/cmake/options.cmake`)
   :header-rows: 1
   :widths: 32 14 54

   * - Option
     - Default
     - Effect
   * - ``CMAKE_BUILD_TYPE``
     - ``Release``
     - CMake's own flags (gcc/clang: Release ``-O3 -DNDEBUG``, Debug ``-g``); ``Debug`` keeps
       ``PHARE_DEBUG_DO`` checks and SAMRAI assertions. Only ``coverage`` adds ``-O0 -g3``
       (:src:`res/cmake/coverage.cmake#L9`).
   * - ``devMode``
     - OFF
     - Warnings become errors (:src:`res/cmake/def.cmake#L74-L80`), errors become exceptions
       (:doc:`hpc_internals`), and each configure runs ``git pull`` in the already cloned
       subprojects except googletest (:src:`res/cmake/dep.cmake#L3-L15`).
   * - ``test`` / ``testMPI``
     - ON / OFF
     - With ``testMPI``: C++ tests under ``mpirun`` with ``PHARE_MPI_PROCS`` ranks (default 2,
       :src:`res/cmake/def.cmake#L179-L184`), serial-only tests skipped, ``phare_python3_exec`` runs on 2.
   * - ``PHARE_EXEC_LEVEL_MIN`` / ``_MAX``
     - 1 / 10
     - Which tests exist (below).
   * - ``asan`` / ``ubsan`` / ``bench``
     - OFF
     - Sanitizers; micro-benchmarks.
   * - ``withPhlop`` / ``withCaliper``
     - OFF
     - Timers, see :doc:`/hpc/measuring_performance`.
   * - ``-DPHARE_DIAG_DOUBLES=1`` (in ``CMAKE_CXX_FLAGS``)
     - float
     - Diagnostics in double precision, both formats
       (:src:`src/diagnostic/detail/h5writer.hpp#L38`, :src:`src/diagnostic/detail/vtkh5_type_writer.hpp#L137`).

Tests
-----

.. plot::
   :include-source: False
   :caption: Raise PHARE_EXEC_LEVEL_MAX to register heavier tests; MAX = 100 still leaves out levels 101 and 121.

   from phare_draw import *
   from matplotlib.patches import Rectangle
   fig, (ax,) = strip(1, ["Test levels: a test exists if MIN ≤ level ≤ MAX (default MIN = 1, MAX = 10)"], height=3.0)
   lv = [(1, "unit"), (9, "simulator"), (10, "dispersion"),
         (11, "functional"), (21, "MHD"), (101, "large\nMHD"), (121, "MHD\nconv.")]
   x = {l: i * 2.0 for i, (l, _) in enumerate(lv)}
   ax.plot([-0.6, 12.6], [0, 0], color=C["grey"], lw=1.2)
   for i, (l, name) in enumerate(lv):
       mark(ax, "primal", x[l], 0, color=C["ink"], size=10)
       text(ax, x[l], 0.45, str(l), bold=True)
       text(ax, x[l], -0.45 - 0.6 * (i % 2), name.replace("\n", " "), fs=11)
   for k, (hi, lab, col) in enumerate([(10, "MAX = 10 (default)", C["L1"]), (21, "MAX = 100", C["L2"]), (121, "MAX = 121", C["bad"])]):
       y = 1.0 + 0.45 * k
       arrow(ax, (x[1] - 0.3, y), (x[hi] + 0.3, y), color=col, lw=2, ms=8, style="-")
       text(ax, x[hi] + 0.45, y, lab, color=col, ha="left", fs=11)
   ax.plot([-0.6, 15.5], [-1.4, 2.2], alpha=0)
   finish(fig, [ax])

* Level filter: ``phare_exec`` / ``phare_mpi_python3_exec`` (:src:`res/cmake/def.cmake#L318-L343`);
  plain unit tests count as level 1.
* Labels (only two at this commit): ``physics`` (dispersion) and ``nightly`` (3D particle
  initialization), set with ``phare_test_labels``, Python tests only (:src:`res/cmake/def.cmake#L352-L360`).
* ``PHARE_DRY_RUN=1``: scripts set up but ``initialize``, ``advance`` and ``run`` do nothing
  (:src:`pyphare/pyphare/simulator/simulator.py#L136-L137`); CI uses it on the heavy levels
  (:src:`res/ci/jobs/pull_request/dry_run_functional_tests.sh#L10-L12`).
* Add a test in the folder's ``CMakeLists.txt``: ``add_phare_test`` (C++),
  ``phare_python3_exec(level name script.py dir)``, ``phare_mpi_python3_exec(level N name script.py dir)``.

.. code-block:: bash

   ctest -j 8 --output-on-failure        # everything registered
   ctest -R restarts                     # names matching
   ctest -L physics                      # one label
   ctest -LE "nightly|physics"           # all but these labels

Continuous integration
----------------------

* **Gate** (``gate.yml``), diff only, no build: skipped tests, loosened tolerances, changed input
  keywords, clang-format, ruff, removed tests, PR description, cyclic imports, clang-tidy, copy-paste,
  workflows, dependencies. Scripts come from ``master`` (:src:`.github/workflows/gate.yml#L1-L14`).
* **Builds**: GitHub Actions with gcc and clang (``cmake_ubuntu.yml``), coverage
  (``ctest -LE "nightly|physics"``), and the project's TeamCity server (configuration not in the repository).
* **Review threads**: a check fails while a thread started by a person is unresolved (``review-threads.yml``).

Conventions
-----------

* An issue first, one issue per PR, small PRs, tests with the code: ``CONTRIBUTING.md`` and the PR template
  (Issue / What this implements / Design / Tests).
* C++: ``.clang-format`` (version pinned in ``pyproject.toml``). Python: lint rules in ``ruff.toml`` (ruff
  ignores the ``[tool.ruff]`` table of ``pyproject.toml`` while ``ruff.toml`` exists); dependencies in
  ``pyproject.toml``; no top-level ``h5py``, ``mpi4py``, ``scipy.optimize`` imports (``doc/conventions.md``).

.. admonition:: In preparation
   :class: in-preparation

   Typical run time of each test level on a workstation.

Go deeper
---------

* :doc:`how_to_add` — recipes for common changes. · :doc:`cpp_api`

Related
-------

* :doc:`/getting_started/build` — user build instructions. · :doc:`/hpc/measuring_performance`
