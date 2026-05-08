## Performance analyses

### via perf

See: https://man7.org/linux/man-pages/man1/perf-record.1.html

```shell
# build with the following
CMAKE_CXX_FLAGS="-g3 -O3 -march=native -mtune=native -fno-omit-frame-pointer"
CMAKE_BUILD_TYPE="RelWithDebInfo"

# build...

perf record -Tga -F 1000 ./build/src/phare/phare-exe path/to/python/script.py
ls -lah perf.data # check file is new/exists
perf script report flamegraph
```

### via perf in vtune

installation see
https://www.intel.com/content/www/us/en/docs/vtune-profiler/installation-guide/2023-0/package-managers.html


```shell
# see above, but rename perf.data to data.perf (vtune only reads .perf files)
```

### via scalasca / scorep

```shell
# install
sudo dnf install scorep-openmpi scalasca-openmpi libunwind-devel binutils-devel elfutils-devel
```

```shell
# usage
CMAKE_CXX_FLAGS="-g3 -O3 -march=native -mtune=native -fno-omit-frame-pointer"
CMAKE_BUILD_TYPE="RelWithDebInfo"
export CXX=scorep-mpicxx CC=scorep-mpicc FC=scorep-gfortran

# build...

# run
# see https://vampir.eu/public/files/pdf/spcheatsheet_a4.pdf
export SCOREP_EXPERIMENT_DIRECTORY=scorep_run_trace
export SCOREP_TOTAL_MEMORY=500M # default is 16MB which will fail
scalasca -analyze mpirun -n 8 ./build/src/phare/phare-exe path/to/python/script.py

# view with GUI
scalasca -examine scorep_run_trace

# or write scorp.score file if no GUI available
scalasca -examine -s scorep_run_trace
```

### via phlop

PHARE exposes a number of logging functions which can be configured for doing stack based scope timings

This can be used as follows:

1. CMake configuration: <br />
        `-DwithPhlop=ON -DCMAKE_CXX_FLAGS="-DPHARE_LOG_LEVEL=1"`

    PHARE_LOG_LEVEL supports the following values
    - 0 = OFF (default)
    - 1 = LIGHT
    - 2 = MEDIUM
    - 3 = HEAVY

    Each C++ scope is configured with a value of 1 to 3, to show its "cost". <br />
    e.g. `PHARE_LOG_SCOPE(1, "Simulator::initialize");` <br />

    A value of 1, should try to be such that there is an equal number of calls on all ranks for that scope. <br />
    2 and 3 are provided for general use and 2 should be called less than 3.  <br />
    A default value can be configured by editing `src/core/logger.hpp`
    e.g. `#define PHARE_LOG_LEVEL 1`

2. Runtime

    Use environment variable to activate it `PHARE_SCOPE_TIMING=1` <br />
    0 = OFF (default) <br />
    1 = ON

3. Analysis

    Finally, when your execution is complete, there should be a timings file per rank.  <br />
    These files may be interrogated with some provided python scripts to display the data

#### example

```shell
PHARE_SCOPE_TIMING=1 mpirun -n 3 python3 tests/functional/harris/harris_3d.py

find .phare/
.phare/
.phare/timings
.phare/timings/rank.1.txt
.phare/timings/rank.2.txt
.phare/timings/rank.0.txt

export PYTHONPATH=${PWD}:${PWD}/build:${PWD}/pyphare:${PWD}/subprojects/phlop
python3 tools/python3/phloping.py print_scope_timings -f .phare/timings/rank.0.txt
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
```


The `loss` value indicates the percentage of time not accounted for by nested stack scope times. <br />
It is recommended to minimize this value.

```bash
python3 tools/python3/phloping.py print_variance_across -f ".phare/timings/rank.*.txt"
 Simulator::initialize 20843 21889
  HybridLevelInitializer::initialize_level 165859 235571
  HybridLevelInitializer::initialize_level 71650 4627296
  HybridLevelInitializer::initialize_level 46113 8926704
 DiagnosticsManager::dump 8844 1424330
 Simulator::advance 1414933 1414505
  HybridLevelInitializer::initialize_level 3612120 41361203
  HybridLevelInitializer::initialize_level 659645 42451743
 DiagnosticsManager::dump 3882 1995649
 Simulator::advance 1984217 1983348
  HybridLevelInitializer::initialize_level 1013733 5362220
  HybridLevelInitializer::initialize_level 1020942 8892336
  HybridLevelInitializer::initialize_level 77815 18638902
 DiagnosticsManager::dump 2218 2054245
```

Here each line prints the standard deviation for the scope start time, and run time for all ranks, in nanoseconds <br />
This feature only works if there is a constant number of scope calls across all ranks. <br />
i.e. there are no patch based loop stack logs which would likely vary across ranks.
