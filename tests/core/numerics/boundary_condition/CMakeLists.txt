cmake_minimum_required (VERSION 3.20.1)

project(test-boundary_condition)

set(SOURCES test_main.cpp)

add_executable(${PROJECT_NAME} ${SOURCES})

target_include_directories(${PROJECT_NAME} PRIVATE
  ${GTEST_INCLUDE_DIRS}
  )

target_link_libraries(${PROJECT_NAME} PRIVATE
  phare_core
  cppdict
  ${GTEST_LIBS})

add_no_mpi_phare_test(${PROJECT_NAME} ${CMAKE_CURRENT_BINARY_DIR})


foreach(suffix none neumann dirichlet symmetric antisymmetric
               div_free_transverse_neumann
               div_free_transverse_dirichlet
               dirichlet_b_energy_subbc
               total_energy_from_pressure
               inflow_compose)
  set(target test-field-bc-${suffix})
  add_executable(${target} test_field_boundary_conditions_${suffix}.cpp)
  target_include_directories(${target} PRIVATE ${GTEST_INCLUDE_DIRS} ${PHARE_PROJECT_DIR})
  target_link_libraries(${target} PRIVATE phare_core cppdict ${GTEST_LIBS})
  add_no_mpi_phare_test(${target} ${CMAKE_CURRENT_BINARY_DIR})
endforeach()


