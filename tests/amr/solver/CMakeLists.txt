cmake_minimum_required (VERSION 3.20.1)

project(test-solver-ppc)

add_executable(${PROJECT_NAME} test_solver_ppc.cpp)

target_include_directories(${PROJECT_NAME} PRIVATE
  ${GTEST_INCLUDE_DIRS}
  )

target_link_libraries(${PROJECT_NAME} PRIVATE
  phare_amr
  phare_simulator
  ${GTEST_LIBS})

add_no_mpi_phare_test(${PROJECT_NAME} ${CMAKE_CURRENT_BINARY_DIR})
