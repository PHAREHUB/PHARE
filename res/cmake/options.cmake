# PHARE CMake options file
#  If adding new option, consider adding it to the function "print_phare_options"


# 1. Highfive
option(HighFive "Build with highfive usage" ON)
# Downloads Highfive if missing from: https://github.com/BlueBrain/HighFive

# 2. test
option(test "Build test with google test" ON)
# Configures targets used with the fuction "add_phare_test" to compile and be tested.
# Testing is run serially and locally without MPI

# 3. testMPI
option(testMPI "Run tests in parallel with mpiriun" OFF)
# Overrides test procedure from -Dtest=ON
# Configures all tests to run via "mpirun -n $PHARE_MPI_PROCS"
# PHARE_MPI_PROCS is a cmake argument which if unset default is "2"

# -Dcoverage=ON
option(coverage "Generate coverage" OFF)
# Enables coverage and generation of coverage from tests via gcovr


# -DdevMode=ON
option(devMode "Build with Werror/etc" OFF)
# Enables stricter compiler flags

# -Dcppcheck=ON
option(cppcheck "Enable cppcheck xml report" OFF)
#

# -Dasan=ON
option(asan "build with asan support" OFF)
#

# -Dubsan=ON
option(ubsan "build with ubsan support" OFF)
#

# -DforceSerialTest=ON
option(forceSerialTest "force the test to be serial under MPI mode" OFF)
# This flag gives the option to still build/test targets which are defined via the functions:
#    add_no_mpi_phare_test
#    add_no_mpi_python3_test


# -DforceGetPybind=ON
option(forceGetPybind "force retrieval of pybind from github" OFF)
# Useful if you want to avoid using the system wide installed version of pybind


#  -DwithCcache=ON # doesn't include dependencies
option(withCcache "build with ccache if found and devMode=ON" ON)


# -DwithIPO=OFF
option(withIPO "Use IPO/LTO if system supported" OFF)
# IPO can cause linking to fail so default off
# https://cmake.org/cmake/help/latest/module/CheckIPOSupported.html

# -DwithCaliper=OFF
option(withCaliper "Use LLNL Caliper" OFF)


# -DwithPhlop=OFF
option(withPhlop "Use phlop" OFF)


# -DlowResourceTests=ON
option(lowResourceTests "Disable heavy tests for CI (2d/3d/etc" OFF)

# -DtestDuringBuild=ON enabled if devMode=ON, disabled if asan=ON (needs LD_PRELOAD)
option(testDuringBuild "Runs C++ unit tests after they are built" OFF)


# -DPGO_GEN=OFF profile guided optimization generate
option(PGO_GEN "profile guided optimization generate" OFF)
# -DPGO_USE=OFF
option(PGO_USE "profile guided optimization use" OFF)


# Controlling the activation of tests
if (NOT DEFINED PHARE_EXEC_LEVEL_MIN)
  set(PHARE_EXEC_LEVEL_MIN 1)
endif()
if (NOT DEFINED PHARE_EXEC_LEVEL_MAX)
  set(PHARE_EXEC_LEVEL_MAX 10)
endif()

# -Dbench=OFF
option(bench "Compile PHARE Benchmarks" OFF)


# -Dphare_configurator=OFF
option(phare_configurator "Best guess setup for compile" OFF)


# print options
function(print_phare_options)

  message("PHARE CMAKE OPTION LIST ")
  message("Build with highfive usage                   : " ${HighFive})
  message("Build with strict flags e.g. Werror         : " ${devMode})
  message("Build test with google test                 : " ${test})
  message("Build bench with google benchmark           : " ${bench})
  message("Run test with MPI                           : " ${testMPI})
  message("Generate coverage                           : " ${coverage})
  message("Enable cppcheck xml report                  : " ${cppcheck})
  message("force the test to be serial under MPI mode  : " ${forceSerialTest})
  message("Enable IPO/LTO if system supported          : " ${withIPO})
  message("build with ubsan support                    : " ${ubsan})
  message("build with asan support                     : " ${asan})
  message("build with ccache (if found) in devMode     : " ${withCcache})
  message("build with LLNL Caliper                     : " ${withCaliper})
  message("profile guided optimization generate        : " ${PGO_GEN})
  message("profile guided optimization use             : " ${PGO_USE})

  message("MPI_LIBRARY_PATH                            : " ${MPI_LIBRARY_PATH})

  if(${devMode})
    message("PHARE_EXEC_LEVEL_MIN                        : " ${PHARE_EXEC_LEVEL_MIN})
    message("PHARE_EXEC_LEVEL_MAX                        : " ${PHARE_EXEC_LEVEL_MAX})
    message("ccache found                                : " ${PHARE_WITH_CCACHE})
    message("PHARE_BASE_LIBS                             : " ${PHARE_BASE_LIBS})
    message("Exclude heavier tests                       : " ${lowResourceTests})
    message("PHARE_MPI_PROCS                             : " ${PHARE_MPI_PROCS})
    message("Selected Compiler                           : " ${CMAKE_CXX_COMPILER_ID})
  endif()

endfunction(print_phare_options)
