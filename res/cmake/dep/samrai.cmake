
find_package(MPI REQUIRED COMPONENTS C)
# some linkers (like mold) don't use default library paths
get_filename_component(MPI_LIBRARY_PATH ${MPI_LIBRARY} DIRECTORY)


if(DEFINED SAMRAI_ROOT)
  find_package(SAMRAI PATHS ${SAMRAI_ROOT} REQUIRED)
else()

  if(NOT DEFINED PHARE_SAMRAI_VERSION)
    SET(PHARE_SAMRAI_VERSION "develop")
  endif()

  set(SAMRAI_SRCDIR ${CMAKE_CURRENT_SOURCE_DIR}/subprojects/samrai)
  set(SAMRAI_BIN ${CMAKE_CURRENT_BINARY_DIR}/subprojects/samrai)

  phare_github_get_or_update(SAMRAI ${SAMRAI_SRCDIR} LLNL/SAMRAI ${PHARE_SAMRAI_VERSION})

  if (SAMRAI_DEBUG)
    set (DEBUG_INITIALIZE_UNDEFINED Off CACHE BOOL "YES!" FORCE)
    set (DEBUG_CHECK_ASSERTIONS Off CACHE BOOL "YES!" FORCE)
    set (DEBUG_CHECK_DIM_ASSERTIONS Off CACHE BOOL "YES!" FORCE)
  else()
    set (DEBUG_INITIALIZE_UNDEFINED Off CACHE BOOL "NO!" FORCE)
    set (DEBUG_CHECK_ASSERTIONS Off CACHE BOOL "NO!" FORCE)
    set (DEBUG_CHECK_DIM_ASSERTIONS Off CACHE BOOL "NO!" FORCE)
  endif(SAMRAI_DEBUG)

  if(SAMRAI_BUILD_SHARED_LIBS)
    set(BUILD_SHARED_LIBS ON CACHE BOOL "Make shared libs" FORCE) # default as of 25/JAN/2026 is static libs
  endif(SAMRAI_BUILD_SHARED_LIBS)

  set(ENABLE_TIMERS OFF CACHE BOOL "NO!" FORCE)
  set(ENABLE_CHECK_ASSERTIONS OFF CACHE BOOL "NO!" FORCE)

  option(ENABLE_TESTS "Enable Samrai Test" OFF ) # disable SAMRAI Test so that we can use the googletest pulled after
  option(ENABLE_SAMRAI_TESTS "Enable Samrai Test" OFF ) # disable SAMRAI Test so that we can use the googletest pulled after

  add_subdirectory(${SAMRAI_SRCDIR})
  unset(CMAKE_RUNTIME_OUTPUT_DIRECTORY CACHE) # undoes what samrai does, so ctest can continue to work
  include_directories(${CMAKE_BINARY_DIR}/include) # this is needed to find build-dir/include/SAMRAI/SAMRAI_config.h

endif(DEFINED SAMRAI_ROOT)
