

set (PHARE_HAS_HIGHFIVE "0")
if(HighFive)
  # Setup HDF5 first or Highfive will
  message("HighFive enabled - checking HDF5")

  if(DEFINED HDF5_ROOT)
    find_package(HDF5 PATHS ${HDF5_ROOT} REQUIRED)
  else()
    find_package(HDF5 REQUIRED)
  endif()

  message(STATUS "HDF5_LIBRARIES " ${HDF5_LIBRARIES})
  message(STATUS "HDF5_INCLUDE_DIRS " ${HDF5_INCLUDE_DIRS})
  message(STATUS "HDF5_LIBRARY_PATH " ${HDF5_LIBRARY_PATH})
  include_directories(${HDF5_INCLUDE_DIRS})  # clangd mostly

  if(NOT DEFINED PHARE_HIGHFIVE_VERSION)
    SET(PHARE_HIGHFIVE_VERSION "main")
  endif()

  set (HIGHFIVE_SRC ${CMAKE_CURRENT_SOURCE_DIR}/subprojects/highfive)

  phare_github_get_or_update(HighFive ${HIGHFIVE_SRC} highfive-devs/highfive ${PHARE_HIGHFIVE_VERSION})

  set(HIGHFIVE_UNIT_TESTS OFF) # silence warning
  set(HIGHFIVE_USE_BOOST OFF)
  set(HIGHFIVE_BUILD_DOCS OFF) # conflicts with phare doc target
  set(HIGHFIVE_EXAMPLES OFF)
  include_directories(${HIGHFIVE_SRC}/include)
  add_subdirectory(${HIGHFIVE_SRC})

  if(DEFINED HDF5_ENABLE_PARALLEL AND "${HDF5_ENABLE_PARALLEL}" STREQUAL "ON")
    set (HDF5_IS_PARALLEL TRUE) # this flag is needed if hdf5 is built from source.
  endif()

  if(${HDF5_IS_PARALLEL})
      message("HDF5 PARALLEL detected")
      add_definitions(-DH5_HAVE_PARALLEL)
  else()
      message(WARNING "HDF5 NOT PARALLEL")
  endif()

  set (PHARE_HAS_HIGHFIVE "1")
endif()
