
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

  if (NOT CMAKE_BUILD_TYPE STREQUAL "Release")
    # enable samrai assertions if not in release mode
    set (DEBUG_INITIALIZE_UNDEFINED On)
    set (DEBUG_CHECK_ASSERTIONS On)
    set (DEBUG_CHECK_DIM_ASSERTIONS On)
  endif()

  # CRITICAL: Force SAMRAI to build as shared libraries
  # Static SAMRAI libraries cause multiple instances of static singletons when
  # loaded by different Python modules (cpp_1_2_3, cpp_1_2_4, etc.)
  # This leads to crashes and undefined behavior. See ISSUES.TXT #3
  if(SAMRAI_BUILD_SHARED_LIBS)
    option(BUILD_SHARED_LIBS "Build SAMRAI as shared libraries" ON ) # default as of 15/JUL/2025 is static libs
  endif(SAMRAI_BUILD_SHARED_LIBS)

  option(ENABLE_TESTS "Enable Samrai Test" OFF ) # disable SAMRAI Test so that we can use the googletest pulled after
  option(ENABLE_SAMRAI_TESTS "Enable Samrai Test" OFF ) # disable SAMRAI Test so that we can use the googletest pulled after

  add_subdirectory(${SAMRAI_SRCDIR})
  unset(CMAKE_RUNTIME_OUTPUT_DIRECTORY CACHE) # undoes what samrai does, so ctest can continue to work
  include_directories(${CMAKE_BINARY_DIR}/include) # this is needed to find build-dir/include/SAMRAI/SAMRAI_config.h


  # ============================================================================
  # TEMPORARY WORKAROUND - REMOVE AFTER https://github.com/LLNL/SAMRAI/pull/294
  # These explicit linkages fix missing dependencies in SAMRAI's CMake config
  # ============================================================================
  target_link_libraries(SAMRAI_algs PUBLIC SAMRAI_mesh)
  target_link_libraries(SAMRAI_appu PUBLIC SAMRAI_geom)
  target_link_libraries(SAMRAI_mesh PUBLIC SAMRAI_pdat)
  target_link_libraries(SAMRAI_solv PUBLIC SAMRAI_geom)


endif(DEFINED SAMRAI_ROOT)
