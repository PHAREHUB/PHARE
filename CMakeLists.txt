cmake_minimum_required (VERSION 3.3)

project(PHARE VERSION 0.1)

set(CMAKE_CXX_STANDARD 17)
if (POLICY CMP0074) # hides warning about ${PACKAGE}_ROOT variables
  cmake_policy(SET CMP0074 NEW)
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "" FORCE)
endif()

set(PHARE_PROJECT_DIR ${CMAKE_CURRENT_SOURCE_DIR})

option(HighFive "Build with highfive usage" ON)
option(test "Build test with google test" ON)
option(testMPI "Run tests in parallel with mpiriun" OFF) # -n = PHARE_MPI_PROCS or 2 if unset
option(coverage "Generate coverage" OFF)
option(documentation "Add doxygen target to generate documentation" OFF)
option(cppcheck "Enable cppcheck xml report" OFF)
option(asan "build with asan support" OFF)
option(ubsan "build with ubsan support" OFF)
option(msan "build with msan support" OFF)

find_program(Git git)

include(CheckCXXCompilerFlag)

find_package (PythonInterp  3.0 REQUIRED)

set(SUBPROJECTS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/subprojects)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

#* variables/functions
include("${PHARE_PROJECT_DIR}/res/cmake/def.cmake")
include("${PHARE_PROJECT_DIR}/res/cmake/mkn.kul.cmake")

#* Cppcheck option
include("${PHARE_PROJECT_DIR}/res/cmake/cppcheck.cmake")

#* SAMRAI build option
include("${PHARE_PROJECT_DIR}/res/cmake/samrai.cmake")

#* pybind
include("${PHARE_PROJECT_DIR}/res/cmake/pybind.cmake")

#* HighFive
include("${PHARE_PROJECT_DIR}/res/cmake/highfive.cmake")

include_directories(${PHARE_PROJECT_DIR}/src)

#* coverage and Test option
include("${PHARE_PROJECT_DIR}/res/cmake/coverage.cmake")
include("${PHARE_PROJECT_DIR}/res/cmake/test.cmake")

#*******************************************************************************
#* Build the different libs and executables
#*******************************************************************************

if(NOT ${PHARE_PROJECT_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  file(COPY ${PHARE_PROJECT_DIR}/phare DESTINATION ${CMAKE_BINARY_DIR})
endif()

add_subdirectory(src/core)
add_subdirectory(src/solver)
add_subdirectory(src/diagnostic)
add_subdirectory(src/initializer)
add_subdirectory(src/phare)
add_subdirectory(src/amr)
add_subdirectory(src/simulator)

#* Documentation option
include("${PHARE_PROJECT_DIR}/res/cmake/doc.cmake")

message("PHARE CMAKE OPTION LIST ")
message("Build test with google test                 : " ${test})
message("Run test with MPI                           : " ${testMPI})
message("Generate coverage                           : " ${coverage})
message("Enable cppcheck xml report                  : " ${cppcheck})
message("build with asan support                     : " ${asan})
message("build with ubsan support                    : " ${ubsan})
message("build with msan support                     : " ${msan})
message("Add doxygen target to generate documentation: " ${documentation})
message("Build with highfive usage                   : " ${HighFive})

