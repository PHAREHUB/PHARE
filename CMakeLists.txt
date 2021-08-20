cmake_minimum_required (VERSION 3.9)

project(PHARE VERSION 0.1)
set(PHARE_PROJECT_DIR ${CMAKE_CURRENT_SOURCE_DIR})

include("${PHARE_PROJECT_DIR}/res/cmake/policies.cmake" NO_POLICY_SCOPE)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "" FORCE)
endif()

include("${PHARE_PROJECT_DIR}/res/cmake/options.cmake")
include("${PHARE_PROJECT_DIR}/res/cmake/def.cmake")
include("${PHARE_PROJECT_DIR}/res/cmake/coverage.cmake")
include("${PHARE_PROJECT_DIR}/res/cmake/dep.cmake")

# set standard after dependencies to avoid potential conflicts
set(CMAKE_CXX_STANDARD 17)

# Enables -fPIC  typically required for shared libraries with dependencies. e.g. pybind modules
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

include(CheckCXXCompilerFlag)
find_package(PythonInterp  3.0 REQUIRED)

set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

include_directories(${PHARE_PROJECT_DIR}/src)

include("${PHARE_PROJECT_DIR}/res/cmake/cppcheck.cmake")

set(CMAKE_CXX_STANDARD 17)

#### post deps config
if(PHARE_WITH_CCACHE)
  set(CMAKE_CXX_COMPILER_LAUNCHER ccache)
endif()
####

include("${PHARE_PROJECT_DIR}/res/cmake/test.cmake")
include("${PHARE_PROJECT_DIR}/res/cmake/bench.cmake")

#*******************************************************************************
# PHARE source code
add_subdirectory(src/core)
add_subdirectory(src/initializer)
add_subdirectory(src/amr)
add_subdirectory(src/diagnostic)
add_subdirectory(src/simulator)
add_subdirectory(src/python3)
add_subdirectory(src/phare)
#*******************************************************************************

#* Documentation option
include("${PHARE_PROJECT_DIR}/res/cmake/doc.cmake")

print_phare_options()
