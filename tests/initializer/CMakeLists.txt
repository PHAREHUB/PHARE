cmake_minimum_required (VERSION 3.3)

project(test-initializer)

set(SOURCES test_initializer.cpp)

add_executable(${PROJECT_NAME} ${SOURCES})

target_include_directories(${PROJECT_NAME} PRIVATE
  ${GTEST_INCLUDE_DIRS}
  )


target_link_libraries(${PROJECT_NAME} PRIVATE
  phare_initializer phare_core
  ${GTEST_LIBS})

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/init.py ${CMAKE_CURRENT_BINARY_DIR}/init.py @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/job.py ${CMAKE_CURRENT_BINARY_DIR}/job.py @ONLY)

add_phare_test(${PROJECT_NAME} ${CMAKE_CURRENT_BINARY_DIR})


