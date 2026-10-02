cmake_minimum_required (VERSION 3.20.1)

project(test-vecfield)

set(SOURCES test_main.cpp)

add_executable(${PROJECT_NAME} ${SOURCES})

target_include_directories(${PROJECT_NAME} PRIVATE
  ${GTEST_INCLUDE_DIRS}
  )

target_link_libraries(${PROJECT_NAME} PRIVATE
  phare_core
  ${GTEST_LIBS})

add_no_mpi_phare_test(${PROJECT_NAME} ${CMAKE_CURRENT_BINARY_DIR})



add_executable(test-vecfield-initializer test_vecfield_initializer.cpp)

target_include_directories(test-vecfield-initializer PRIVATE
  ${GTEST_INCLUDE_DIRS}
  )

target_link_libraries(test-vecfield-initializer PRIVATE
  phare_core
  phare_initializer
  ${GTEST_LIBS})

add_no_mpi_phare_test(test-vecfield-initializer ${CMAKE_CURRENT_BINARY_DIR})
