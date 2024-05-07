cmake_minimum_required (VERSION 3.20.1)

project(test-ampere)

set(SOURCES test_main.cpp)

add_executable(${PROJECT_NAME} ${SOURCES})

target_include_directories(${PROJECT_NAME} PRIVATE
  ${GTEST_INCLUDE_DIRS}
  )

target_link_libraries(${PROJECT_NAME} PRIVATE
  phare_core
  ${GTEST_LIBS})

add_custom_command(TARGET test-ampere
                   PRE_BUILD
                   COMMAND  "PYTHONPATH=${PHARE_PYTHONPATH}" ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_ampere.py ${CMAKE_CURRENT_BINARY_DIR})

add_no_mpi_phare_test(${PROJECT_NAME} ${CMAKE_CURRENT_BINARY_DIR})
