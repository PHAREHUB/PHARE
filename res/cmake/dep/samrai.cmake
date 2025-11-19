
find_package(MPI REQUIRED COMPONENTS C)
# some linkers (like mold) don't use default library paths
get_filename_component(MPI_LIBRARY_PATH ${MPI_LIBRARY} DIRECTORY)


find_package(SAMRAI CONFIG QUIET)
if (NOT SAMRAI_FOUND)
    message("SAMRAI NOT FOUND")
  if(DEFINED SAMRAI_ROOT)
    find_package(SAMRAI PATHS ${SAMRAI_ROOT} REQUIRED)
  else()

    if(NOT DEFINED PHARE_SAMRAI_VERSION)
      SET(PHARE_SAMRAI_VERSION "develop")
      SET(PHARE_SAMRAI_VERSION "feature/srcmask") # TORM
    endif()

    set(SAMRAI_SRCDIR ${CMAKE_CURRENT_SOURCE_DIR}/subprojects/samrai)
    set(SAMRAI_BIN ${CMAKE_CURRENT_BINARY_DIR}/subprojects/samrai)

    # phare_github_get_or_update(SAMRAI ${SAMRAI_SRCDIR} LLNL/SAMRAI ${PHARE_SAMRAI_VERSION}) # uncomment
    phare_github_get_or_update(SAMRAI ${SAMRAI_SRCDIR} nicolasaunai/SAMRAI ${PHARE_SAMRAI_VERSION}) # TORM

    if (NOT CMAKE_BUILD_TYPE STREQUAL "Release")
      # enable samrai assertions if not in release mode
      set (DEBUG_INITIALIZE_UNDEFINED On)
      set (DEBUG_CHECK_ASSERTIONS On)
      set (DEBUG_CHECK_DIM_ASSERTIONS On)
    endif()

    option(ENABLE_TESTS "Enable Samrai Test" OFF ) # disable SAMRAI Test so that we can use the googletest pulled after
    option(ENABLE_SAMRAI_TESTS "Enable Samrai Test" OFF ) # disable SAMRAI Test so that we can use the googletest pulled after

    add_subdirectory(${SAMRAI_SRCDIR})
    unset(CMAKE_RUNTIME_OUTPUT_DIRECTORY CACHE) # undoes what samrai does, so ctest can continue to work
    include_directories(${CMAKE_BINARY_DIR}/include) # this is needed to find build-dir/include/SAMRAI/SAMRAI_config.h
  endif()
else()
  include_directories(${SAMRAI_INCLUDE_DIRS}) # this is needed to find build-dir/include/SAMRAI/SAMRAI_config.h
  message("SAMRAI HAS BEEN FOUND")
  message(${SAMRAI_INCLUDE_DIRS})
endif()

