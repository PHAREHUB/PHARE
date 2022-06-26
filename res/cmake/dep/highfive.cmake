

set (PHARE_HAS_HIGHFIVE "0")
if(HighFive)

  set (HIGHFIVE_SRC ${CMAKE_CURRENT_SOURCE_DIR}/subprojects/highfive)

  if (NOT EXISTS ${HIGHFIVE_SRC})
    execute_process(
      COMMAND ${Git} clone https://github.com/BlueBrain/HighFive ${HIGHFIVE_SRC} -b master --depth 1
    )

    # Remove symlinks that cause infinite looping
    #  https://github.com/BlueBrain/HighFive/issues/338
    file(REMOVE_RECURSE ${HIGHFIVE_SRC}/tests)

  else()
    if(devMode)
      message("downloading latest HighFive updates")
      execute_process(COMMAND ${Git} checkout tests WORKING_DIRECTORY ${HIGHFIVE_SRC}) # undo tests delete
      execute_process(COMMAND ${Git} pull origin master WORKING_DIRECTORY ${HIGHFIVE_SRC})
      file(REMOVE_RECURSE ${HIGHFIVE_SRC}/tests) # redelete symlinks
    endif(devMode)
  endif()

  include_directories(
    ${HIGHFIVE_SRC}/include
    ${CMAKE_BINARY_DIR}/subprojects/highfive/include # configured include for version info
  )
  set(HIGHFIVE_UNIT_TESTS OFF) # silence warning
  set(HIGHFIVE_USE_BOOST OFF)
  set(HIGHFIVE_BUILD_DOCS OFF) # conflicts with phare doc target
  set(HIGHFIVE_EXAMPLES OFF)
  add_subdirectory(${HIGHFIVE_SRC})

  message("HighFive enabled - checking HDF5")

  if(DEFINED HDF5_ROOT)
    find_package(HDF5 PATHS ${HDF5_ROOT} REQUIRED)
  else()
    find_package(HDF5 REQUIRED)
  endif()

  if(DEFINED HDF5_ENABLE_PARALLEL AND "${HDF5_ENABLE_PARALLEL}" STREQUAL "ON")
    set (HDF5_IS_PARALLEL TRUE) # this flag is needed if hdf5 is built from source.
  endif()

  if(${HDF5_IS_PARALLEL})
      message("HDF5 PARALLEL detected")
      add_definitions(-DH5_HAVE_PARALLEL)
  else()
      message(WARNING "HDF5 NOT PARALLEL")
  endif()

  message(STATUS "HDF5_LIBRARIES " ${HDF5_LIBRARIES})
  message(STATUS "HDF5_INCLUDE_DIRS " ${HDF5_INCLUDE_DIRS})
  message(STATUS "HDF5_LIBRARY_PATH " ${HDF5_LIBRARY_PATH})
  set (PHARE_HAS_HIGHFIVE "1")
endif()
