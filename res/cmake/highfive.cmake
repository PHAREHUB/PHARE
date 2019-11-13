
if(HighFive)
  include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/subprojects/highfive/include
  )

  if(DEFINED HDF5_ROOT)
    find_package(HDF5 PATHS ${HDF5_ROOT} REQUIRED)
  else()
    find_package(HDF5 REQUIRED)
  endif()

  message("HighFive enabled - checking HDF5")
  if(${HDF5_IS_PARALLEL})
      message("HDF5 PARALLEL detected")
      add_definitions(-DH5_HAVE_PARALLEL)
  else()
      message(WARNING "HDF5 NOT PARALLEL")
  endif()

  message(STATUS "HDF5_LIBRARIES " ${HDF5_LIBRARIES})
  message(STATUS "HDF5_INCLUDE_DIRS " ${HDF5_INCLUDE_DIRS})
  message(STATUS "HDF5_LIBRARY_PATH " ${HDF5_LIBRARY_PATH})
endif()
