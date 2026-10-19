cmake_minimum_required (VERSION 3.20.1)

project(test-particles-data-refine)


   

set(SOURCES_CPP
  test_particle_data_refine.cpp
  test_split.cpp
   )

add_executable(${PROJECT_NAME} ${SOURCES_CPP})


target_include_directories(${PROJECT_NAME} PRIVATE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  ${GTEST_INCLUDE_DIRS}
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  )

target_link_libraries(${PROJECT_NAME} PRIVATE
  phare_amr
  ${GTEST_LIBS})


add_no_mpi_phare_test(${PROJECT_NAME} ${CMAKE_CURRENT_BINARY_DIR})




add_executable(test-split-collisions test_split_collisions.cpp)
target_include_directories(test-split-collisions PRIVATE ${GTEST_INCLUDE_DIRS})
target_link_libraries(test-split-collisions PRIVATE phare_amr ${GTEST_LIBS})
add_no_mpi_phare_test(test-split-collisions ${CMAKE_CURRENT_BINARY_DIR})
