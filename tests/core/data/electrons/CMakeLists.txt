cmake_minimum_required (VERSION 3.20.1)

project(test-electrons)

set(SOURCES test_electrons.cpp)

add_executable(${PROJECT_NAME} ${SOURCES})

target_include_directories(${PROJECT_NAME} PRIVATE
  ${GTEST_INCLUDE_DIRS}
  )

target_link_libraries(${PROJECT_NAME} PRIVATE
  phare_core
  ${GTEST_LIBS})

add_phare_test(${PROJECT_NAME} ${CMAKE_CURRENT_BINARY_DIR})


