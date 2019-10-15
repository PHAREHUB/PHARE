cmake_minimum_required (VERSION 3.3)

project(test-initializer)

set(SOURCES test_initializer.cpp)

add_executable(${PROJECT_NAME} ${SOURCES})

target_include_directories(${PROJECT_NAME} PRIVATE
  $<BUILD_INTERFACE:${gtest_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${gmock_SOURCE_DIR}/include>
  )

target_compile_options(${PROJECT_NAME} PRIVATE ${PHARE_WERROR_FLAGS})
target_link_libraries(${PROJECT_NAME} PRIVATE
  phare_initializer phare_core
  gtest
  gmock)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/setpythonpath.py.in ${CMAKE_CURRENT_BINARY_DIR}/setpythonpath.py @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/init.py ${CMAKE_CURRENT_BINARY_DIR}/init.py @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/job.py ${CMAKE_CURRENT_BINARY_DIR}/job.py @ONLY)

add_phare_test(${PROJECT_NAME} ${CMAKE_CURRENT_BINARY_DIR})

include(${PHARE_PROJECT_DIR}/sanitizer.cmake)
