cmake_minimum_required (VERSION 3.9)

project(test-interpolator)

set(SOURCES test_main.cpp)

add_executable(${PROJECT_NAME} ${SOURCES})

target_include_directories(${PROJECT_NAME} PRIVATE
  ${GTEST_INCLUDE_DIRS}
  )

target_link_libraries(${PROJECT_NAME} PRIVATE
  phare_core
  ${GTEST_LIBS})

add_no_mpi_phare_test(${PROJECT_NAME} ${CMAKE_CURRENT_BINARY_DIR})
add_custom_command(TARGET test-interpolator
                   POST_BUILD
                   COMMAND "PYTHONPATH=${PHARE_PYTHONPATH}" ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/interpolator_test.py ${CMAKE_CURRENT_BINARY_DIR})


if(NOT ${PHARE_PROJECT_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  file(GLOB PYFILES "*.py")
  file(COPY ${PYFILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif()

