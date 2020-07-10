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

# -Dtest=ON
option(coverage "Generate coverage" OFF)
# Enables coverage and generation of coverage from tests via gcovr

# -Ddocumentation=ON
option(documentation "Add doxygen target to generate documentation" OFF)
#

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



# print options
function(print_phare_options)

  message("PHARE CMAKE OPTION LIST ")
  message("Build with highfive usage                   : " ${HighFive})
  message("Build with strict flags e.g. Werror         : " ${devMode})
  message("Build test with google test                 : " ${test})
  message("Run test with MPI                           : " ${testMPI})
  message("Generate coverage                           : " ${coverage})
  message("Enable cppcheck xml report                  : " ${cppcheck})
  message("Add doxygen target to generate documentation: " ${documentation})
  message("force the test to be serial under MPI mode  : " ${forceSerialTest})
  message("force retrieval of pybind from github       : " ${forceGetPybind})
  message("build with asan support                     : " ${asan})
  message("build with ubsan support                    : " ${ubsan})

endfunction(print_phare_options)


