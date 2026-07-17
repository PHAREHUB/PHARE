cmake_minimum_required (VERSION 3.20.1) # released April 8, 2021 - https://www.kitware.com/cmake-3-20-1-available-for-download/

project(PHARE VERSION 0.1 LANGUAGES CXX C)
set(PHARE_PROJECT_DIR ${CMAKE_CURRENT_SOURCE_DIR})

include("${PHARE_PROJECT_DIR}/res/cmake/policies.cmake" NO_POLICY_SCOPE)
include(CheckCXXCompilerFlag)

# Release mode default
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()

find_package(Python 3.10 COMPONENTS Interpreter Development REQUIRED)

if(EXISTS "${PHARE_PROJECT_DIR}/release.ver.txt") # delete file to use newer deps
  include("${PHARE_PROJECT_DIR}/res/cmake/release.cmake")
endif()

include("${PHARE_PROJECT_DIR}/res/cmake/options.cmake")
include("${PHARE_PROJECT_DIR}/res/cmake/def.cmake")
include("${PHARE_PROJECT_DIR}/res/cmake/coverage.cmake")
include("${PHARE_PROJECT_DIR}/res/cmake/dep.cmake")

# Enables -fPIC  typically required for shared libraries with dependencies. e.g. pybind modules
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

include_directories(${PHARE_PROJECT_DIR}/src)

include("${PHARE_PROJECT_DIR}/res/cmake/cppcheck.cmake")

# set standard after dependencies to avoid potential conflicts
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#### post deps config
if(PHARE_WITH_CCACHE)
  set(CMAKE_CXX_COMPILER_LAUNCHER ccache)
endif()
####

include("${PHARE_PROJECT_DIR}/res/cmake/test.cmake")

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

include("${PHARE_PROJECT_DIR}/res/cmake/bench.cmake")

print_phare_options()
