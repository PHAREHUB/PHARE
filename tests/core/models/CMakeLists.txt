cmake_minimum_required (VERSION 3.20.1)

project(test-external-field-updater)

set(TESTS test_external_field_updater_dipole
          test_external_field_updater_factory
          test_external_field_updater_user_defined)

foreach(test ${TESTS})
  string(
    REPLACE 
    "test_external_field_updater_"
    "test-external-field-updater-"
    target 
    ${test}
  )

  add_executable(${target} ${test}.cpp)

  target_include_directories(${target} PRIVATE
    ${GTEST_INCLUDE_DIRS}
    )

  target_link_libraries(${target} PRIVATE
    phare_core
    ${GTEST_LIBS})

  add_no_mpi_phare_test(${target} ${CMAKE_CURRENT_BINARY_DIR})
endforeach()
