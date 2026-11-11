cmake_minimum_required (VERSION 3.20.1)

project(test-updater)

set(SOURCES test_updater.cpp)

add_executable(${PROJECT_NAME} ${SOURCES})

target_include_directories(${PROJECT_NAME} PRIVATE
  ${GTEST_INCLUDE_DIRS}
  )

target_link_libraries(${PROJECT_NAME} PRIVATE
  phare_core
  phare_simulator
  ${GTEST_LIBS})

add_no_mpi_phare_test(${PROJECT_NAME} ${CMAKE_CURRENT_BINARY_DIR})



# the updater on every patch of a SAMRAI level, needs SAMRAI for the same-level neighbours
add_executable(test-updater-on-level test_updater_on_level.cpp)

target_include_directories(test-updater-on-level PRIVATE
  ${GTEST_INCLUDE_DIRS}
  )

target_link_libraries(test-updater-on-level PRIVATE
  phare_core
  phare_amr
  phare_simulator
  ${GTEST_LIBS})

add_no_mpi_phare_test(test-updater-on-level ${CMAKE_CURRENT_BINARY_DIR})
