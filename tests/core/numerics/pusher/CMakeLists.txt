cmake_minimum_required (VERSION 3.20.1)

project(test-pusher)

set(SOURCES test_pusher.cpp)

add_executable(${PROJECT_NAME} ${SOURCES})

target_include_directories(${PROJECT_NAME} PRIVATE
  ${GTEST_INCLUDE_DIRS}
  )

target_link_libraries(${PROJECT_NAME} PRIVATE
  phare_core
  phare_initializer
  ${GTEST_LIBS})

add_custom_command(TARGET test-pusher
                   PRE_BUILD
                   COMMAND "PYTHONPATH=${PHARE_PYTHONPATH}" ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_pusher.py ${CMAKE_CURRENT_BINARY_DIR})

add_no_mpi_phare_test(${PROJECT_NAME} ${CMAKE_CURRENT_BINARY_DIR})


add_executable(test-boris-move-error test_boris_move_error.cpp)

target_include_directories(test-boris-move-error PRIVATE
  ${GTEST_INCLUDE_DIRS}
  )

target_link_libraries(test-boris-move-error PRIVATE
  phare_core
  ${PHARE_BASE_LIBS}
  ${GTEST_LIBS})

add_no_mpi_phare_test(test-boris-move-error ${CMAKE_CURRENT_BINARY_DIR})
