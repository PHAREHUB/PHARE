
if(HighFive)

  if (NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/subprojects/highfive)
    execute_process(
      COMMAND ${Git} clone https://github.com/BlueBrain/HighFive ${CMAKE_CURRENT_SOURCE_DIR}/subprojects/highfive -b master --depth 1
    )
  endif()

  include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/subprojects/highfive/include
  )

  message("HighFive enabled - checking HDF5")

  if(DEFINED HDF5_ROOT)
    find_package(HDF5 PATHS ${HDF5_ROOT} REQUIRED)
  else()
    find_package(HDF5 REQUIRED)
  endif()

  if(DEFINED HDF5_ENABLE_PARALLEL AND "${HDF5_ENABLE_PARALLEL}" STREQUAL "ON")
    set (HDF5_IS_PARALLEL TRUE)

    # This is here for HDF5 versions > 1.12
    #  This should be able to be removed when this is merged to master
    #   https://github.com/BlueBrain/HighFive/pull/311
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DH5Oget_info_vers=1 -DH5O_info_t_vers=1")
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
endif()
