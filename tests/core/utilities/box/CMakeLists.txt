cmake_minimum_required (VERSION 3.20.1)

project(test-box)

set(SOURCES test_box.cpp)

add_executable(${PROJECT_NAME} ${SOURCES})

target_include_directories(${PROJECT_NAME} PRIVATE
  ${GTEST_INCLUDE_DIRS}
  )

target_link_libraries(${PROJECT_NAME} PRIVATE
  phare_core
  ${GTEST_LIBS})

add_no_mpi_phare_test(${PROJECT_NAME} ${CMAKE_CURRENT_BINARY_DIR})



add_executable(test-box-span test_box_span.cpp)
target_include_directories(test-box-span PRIVATE ${GTEST_INCLUDE_DIRS})
target_link_libraries(test-box-span PRIVATE phare_core ${GTEST_LIBS})
add_no_mpi_phare_test(test-box-span ${CMAKE_CURRENT_BINARY_DIR})
