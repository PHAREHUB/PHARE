cmake_minimum_required(VERSION 3.20.1)
project(test-ghost-width-calculator)

add_executable(${PROJECT_NAME} test_ghost_width_calculator.cpp)
target_include_directories(${PROJECT_NAME} PRIVATE ${GTEST_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} PRIVATE phare_core ${GTEST_LIBS})
add_no_mpi_phare_test(${PROJECT_NAME} ${CMAKE_CURRENT_BINARY_DIR})
