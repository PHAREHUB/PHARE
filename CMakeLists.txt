cmake_minimum_required (VERSION 3.20.1) # released April 8, 2021 - https://www.kitware.com/cmake-3-20-1-available-for-download/

project(PHARE VERSION 0.1 LANGUAGES CXX C)

# Release mode default
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()

set(CMAKE_POLICY_VERSION_MINIMUM 3.5) # TORM after https://github.com/LLNL/SAMRAI/issues/266

if (POLICY CMP0074) # hides warning about ${PACKAGE}_ROOT variables
  cmake_policy(SET CMP0074 NEW)
endif()
if (POLICY CMP0069) # allow LTO if available, requires cmake 3.9 (released August 2017)
  cmake_policy(SET CMP0069 NEW)
endif()

# Enables -fPIC  typically required for shared libraries with dependencies. e.g. pybind modules
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

include(CheckCXXCompilerFlag)
find_program(Git git)
find_package(Python 3.8 COMPONENTS Interpreter Development REQUIRED)

set(PHARE_PROJECT_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

if(EXISTS "${PHARE_PROJECT_DIR}/release.ver.txt") # delete file to use newer deps
  include("${PHARE_PROJECT_DIR}/res/cmake/release.cmake")
endif()

include_directories(${PHARE_PROJECT_DIR}/src)

include("${PHARE_PROJECT_DIR}/res/cmake/options.cmake")
include("${PHARE_PROJECT_DIR}/res/cmake/def.cmake")
include("${PHARE_PROJECT_DIR}/res/cmake/coverage.cmake")
include("${PHARE_PROJECT_DIR}/res/cmake/dep.cmake")
include("${PHARE_PROJECT_DIR}/res/cmake/cppcheck.cmake")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#### post deps config
if(PHARE_WITH_CCACHE)
  set(CMAKE_CXX_COMPILER_LAUNCHER ccache)
endif()
####

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

include("${PHARE_PROJECT_DIR}/res/cmake/test.cmake")
include("${PHARE_PROJECT_DIR}/res/cmake/bench.cmake")

print_phare_options()
