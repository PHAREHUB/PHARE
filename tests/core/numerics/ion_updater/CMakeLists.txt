cmake_minimum_required (VERSION 3.9)

project(test-updater)


function(_add_ion_updater_test src_name)

  add_executable(${src_name} ${src_name}.cpp)

  target_include_directories(${src_name} PRIVATE
    ${GTEST_INCLUDE_DIRS}
  )

  target_link_directories(${src_name} PRIVATE
    ${HDF5_LIBRARY_PATH}
  )

  target_link_libraries(${src_name} PRIVATE
    phare_core
    phare_simulator
    ${GTEST_LIBS})

  add_no_mpi_phare_test(${src_name} ${CMAKE_CURRENT_BINARY_DIR})

endfunction(_add_ion_updater_test)

_add_ion_updater_test(test_updater)
_add_ion_updater_test(test_lol_updater)
