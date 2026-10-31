

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
  # NOT include_directories()/add_definitions() here -- those are directory-scoped and leak
  # HDF5's (and, since it's MPI-parallel, MPI's) include dirs into every target below this
  # point, including phare_core/phare_mpi which must stay HDF5/MPI-agnostic. HDF5::HDF5 is
  # instead linked PUBLIC on phare_amr only ("SAMRAI uses HDF5" - see amr/CMakeLists.txt),
  # so its include dirs/definitions reach exactly the targets that link phare_amr and no more.

  if(NOT DEFINED PHARE_HIGHFIVE_VERSION)
    SET(PHARE_HIGHFIVE_VERSION "main")
  endif()

  set (HIGHFIVE_SRC ${CMAKE_CURRENT_SOURCE_DIR}/subprojects/highfive)

  phare_github_get_or_update(HighFive ${HIGHFIVE_SRC} highfive-devs/highfive ${PHARE_HIGHFIVE_VERSION})

  set(HIGHFIVE_UNIT_TESTS OFF) # silence warning
  set(HIGHFIVE_USE_BOOST OFF)
  set(HIGHFIVE_BUILD_DOCS OFF) # conflicts with phare doc target
  set(HIGHFIVE_EXAMPLES OFF)
  add_subdirectory(${HIGHFIVE_SRC})

  if(DEFINED HDF5_ENABLE_PARALLEL AND "${HDF5_ENABLE_PARALLEL}" STREQUAL "ON")
    set (HDF5_IS_PARALLEL TRUE) # this flag is needed if hdf5 is built from source.
  endif()

  if(${HDF5_IS_PARALLEL})
      message("HDF5 PARALLEL detected")
  else()
      message(WARNING "HDF5 NOT PARALLEL")
  endif()

  set (PHARE_HAS_HIGHFIVE "1")
endif()
