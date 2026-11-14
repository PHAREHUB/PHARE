=====================
Measuring performance
=====================

**Goal:** find out where the time and memory of your run go, before
changing the number of :term:`ranks <rank>`, the :term:`patch` sizes or the input.

.. contents:: On this page
   :local:
   :depth: 1


PHARE gives you three tools, from the simplest to the most detailed:

.. list-table::
   :header-rows: 1
   :widths: 25 40 35

   * - Tool
     - Tells you
     - Needs
   * - Time per step
     - how long each time step took
     - nothing
   * - Resource monitor
     - memory and CPU use of each rank, sampled in time
     - the ``phlop`` Python package
   * - Scope timers
     - time spent in each part of the code, per rank
     - a build with ``-DwithPhlop=ON`` and a log level

.. mermaid::

   flowchart TD
     RUN["your run"] --> A["time per step<br/>plot_times"]
     RUN --> B["resource monitor<br/>.phare/stats/"]
     RUN --> C["scope timers<br/>.phare/timings/"]
     C --> D["phloping.py"]

Time per step
-------------

``Simulator.run`` prints the wall time of every step and, at the end, the
mean and total:

.. code-block:: text

   t =  0.50000 - 0.00287sec - total 0:00:00.081165
   mean advance time = 0.0008116531372070313
   total advance time = 0:00:00.081165

``Simulator(sim).run(plot_times=True)`` also saves a plot of the time of
each step to ``timestep_times.png`` (rank 0) when the run ends. A step that is suddenly slower often matches
a regrid or a diagnostic dump.

Resource monitor (memory, CPU)
------------------------------

Install the optional package ``phlop`` at the version PHARE is built with
(``python3 -m pip install phlop==0.0.45``),
then ask for monitoring when you run:

.. code-block:: python

   Simulator(sim).run(monitoring=10)   # one sample every 10 seconds

Or, without changing the script, set the environment variable
``PHARE_SIM_MON=1`` (one sample every second). Each rank writes
``.phare/stats/rank.<rank>.yaml``; here written by ``run(monitoring=1)``:

.. code-block:: yaml

   rank: 0
   cli_args:
     interval: 1
   headers:
   - open_file
   - mem_used_mb
   - cpu_usage
   start: 2026-11-13T21:39:20.259048+00:00
   snapshots:
   - v: 3,248,99.9
   - v: 3,248,89.9

Each snapshot gives the number of open files, the memory used by the
process (resident memory, in MB) and its CPU use in percent. See
:doc:`memory` for a script that reads these files.

.. note::

   Without ``phlop`` installed, monitoring silently does nothing.

Scope timers
------------

PHARE's C++ code is instrumented with named timers around its main
functions (``Simulator::advance``, level initialization, diagnostics
dumps...). They cost nothing unless you compile them in.

1. **Build** with phlop and a log level:

   .. code-block:: bash

       -DwithPhlop=ON -DCMAKE_CXX_FLAGS="-DPHARE_LOG_LEVEL=1"

   ``PHARE_LOG_LEVEL`` selects how many timers are compiled in:

   - ``0`` = OFF (default)
   - ``1`` = LIGHT
   - ``2`` = MEDIUM
   - ``3`` = HEAVY

   Each timer in the code is tagged with a level from 1 to 3; a build
   keeps the timers whose level is at most ``PHARE_LOG_LEVEL``. Higher
   levels add finer, more frequently called timers.

2. **Run** with the environment variable
   :ref:`PHARE_SCOPE_TIMING <phare_scope_timing>` set to ``1``. Each rank
   writes ``.phare/timings/rank.<N>.txt``.

3. **Analyse** the files with the provided script.

Example
^^^^^^^

.. code-block:: shell

    PHARE_SCOPE_TIMING=1 mpirun -n 3 python3 tests/functional/harris/harris_3d.py

    export PYTHONPATH=${PWD}:${PWD}/build:${PWD}/pyphare:${PWD}/subprojects/phlop
    python3 tools/python3/phloping.py print_scope_timings -f .phare/timings/rank.0.txt

.. code-block:: text

    100% loss(95.35) Simulator::advance                      32,810.81ms
     1.72% loss(-) HybridLevelInitializer::initialize_level     565.57ms
     1.75% loss(-) HybridLevelInitializer::initialize_level     573.71ms
     1.18% loss(-) HybridLevelInitializer::initialize_level     387.01ms
    100% loss(95.75) Simulator::advance                      31,469.41ms
     2.57% loss(-) HybridLevelInitializer::initialize_level     807.86ms
     1.69% loss(-) HybridLevelInitializer::initialize_level     531.15ms
    100% loss(4.36) Simulator::initialize                     2,183.69ms
     30.60% loss(-) HybridLevelInitializer::initialize_level    668.19ms
     17.15% loss(-) HybridLevelInitializer::initialize_level    374.61ms
     47.89% loss(-) HybridLevelInitializer::initialize_level  1,045.67ms
    100% loss(-) DiagnosticsManager::dump                        77.77ms
    100% loss(-) DiagnosticsManager::dump                        60.24ms
    100% loss(-) DiagnosticsManager::dump                        53.07ms

How to read it:

* A line starting with ``100%`` is a top-level timed scope; the indented
  lines below it are the timed scopes called inside it.
* The percentage is the share of the top-level scope's time.
* ``loss(x)`` is the share (in %) of this scope's time spent *outside* its
  timed sub-scopes, i.e. in code without a timer. ``-`` means no timed
  sub-scope, or a negligible loss.

.. admonition:: In preparation
   :class: in-preparation

   A worked analysis of a production run: comparing the timers across
   ranks to spot load imbalance.

Other profilers
---------------

PHARE is ordinary C++ and MPI code: general tools such as ``perf`` or
Score-P/Scalasca should work on it, as on any such code. Build with ``CMAKE_BUILD_TYPE=RelWithDebInfo``
and add ``-fno-omit-frame-pointer`` to ``CMAKE_CXX_FLAGS`` so call stacks
are readable. PHARE can also send its timers to the Caliper profiling
library instead of phlop (CMake option ``-DwithCaliper=ON``).

Common mistakes
---------------

* **No timing file** — ``PHARE_SCOPE_TIMING=1`` on a build without
  ``-DwithPhlop=ON`` writes nothing, silently — rebuild with it.
* **No timer in the output** — built with ``-DwithPhlop=ON`` but without
  ``-DPHARE_LOG_LEVEL=1`` (or 2, 3): every timer is compiled out — add it
  to ``CMAKE_CXX_FLAGS``.
* **No** ``.phare/stats/`` **files** — ``phlop`` not installed: monitoring
  silently does nothing — install it.

.. seealso:: How the timers are placed in the code: :doc:`/dev/hpc_internals`.
