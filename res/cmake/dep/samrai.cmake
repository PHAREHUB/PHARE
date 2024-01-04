
find_package(MPI REQUIRED COMPONENTS C)
# some linkers (like mold) don't use default library paths
get_filename_component(MPI_LIBRARY_PATH ${MPI_LIBRARY} DIRECTORY)

find_package(SAMRAI CONFIG QUIET)
if (NOT SAMRAI_FOUND)
    message("SAMRAI NOT FOUND")
  if(DEFINED SAMRAI_ROOT)
    find_package(SAMRAI PATHS ${SAMRAI_ROOT} REQUIRED)
  else()
    set(SAMRAI_SRCDIR ${CMAKE_CURRENT_SOURCE_DIR}/subprojects/samrai)
    set(SAMRAI_BIN ${CMAKE_CURRENT_BINARY_DIR}/subprojects/samrai)

    if (NOT EXISTS ${SAMRAI_SRCDIR})
      execute_process(
        COMMAND ${Git} clone https://github.com/LLNL/SAMRAI ${SAMRAI_SRCDIR} -b develop --recursive --depth 10
        )
    else()
      if(devMode)
        message("downloading latest SAMRAI updates")
        execute_process(COMMAND ${Git} pull origin master WORKING_DIRECTORY ${SAMRAI_SRCDIR})
      endif(devMode)
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

