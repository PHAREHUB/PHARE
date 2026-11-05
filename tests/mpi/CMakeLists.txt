cmake_minimum_required (VERSION 3.20.1)

project(test-mpi-utils)

add_executable(${PROJECT_NAME} test_mpi_utils.cpp)

target_include_directories(${PROJECT_NAME} PRIVATE
  ${GTEST_INCLUDE_DIRS}
)

target_link_libraries(${PROJECT_NAME} PRIVATE
  phare_mpi
  phare_core
  MPI::MPI_C
  ${GTEST_LIBS})

add_phare_test(${PROJECT_NAME} ${CMAKE_CURRENT_BINARY_DIR})
