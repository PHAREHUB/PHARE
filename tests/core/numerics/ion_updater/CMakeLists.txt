cmake_minimum_required (VERSION 3.20.1)

project(test-updater)


function(_add_ion_updater_test src_name test_name)

  add_executable(${test_name} ${src_name}.cpp)

  target_include_directories(${test_name} PRIVATE
    ${GTEST_INCLUDE_DIRS}
  )

  target_link_directories(${test_name} PRIVATE
    ${HDF5_LIBRARY_PATH}
  )

  target_link_libraries(${test_name} PRIVATE
    phare_core
    phare_simulator
    pybind11::embed
    ${GTEST_LIBS})

  add_no_mpi_phare_test(${test_name} ${CMAKE_CURRENT_BINARY_DIR})

endfunction(_add_ion_updater_test)

_add_ion_updater_test(test_updater test-updater)
_add_ion_updater_test(test_multi_updater test_multi_updater)
#_add_ion_updater_test(test_multithreaded_updater test_multithreaded_updater)
