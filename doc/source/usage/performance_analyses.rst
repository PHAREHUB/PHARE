=====================
Performance analyses
=====================

via phlop
---------

PHARE exposes a number of logging functions which can be configured for doing
stack based scope timings.

This can be used as follows:

1. CMake configuration

.. code-block:: bash

    -DwithPhlop=ON -DCMAKE_CXX_FLAGS="-DPHARE_LOG_LEVEL=1"

``PHARE_LOG_LEVEL`` supports the following values:

- ``0`` = OFF (default)
- ``1`` = LIGHT
- ``2`` = MEDIUM
- ``3`` = HEAVY

2. Runtime

Use the environment variable :ref:`PHARE_SCOPE_TIMING <phare_scope_timing>`
to activate it:

- ``0`` = OFF (default)
- ``1`` = ON

3. Analysis

Finally, when your execution is complete, there should be a timings file per
rank. These files may be interrogated with some provided python scripts to
display the data.

Example
-------

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
