cmake_minimum_required (VERSION 3.20.1)

project(test-particles-data-stream-pack-unpack)


set(SOURCES_INC
   )

set(SOURCES_CPP
  test_main.cpp
   )

add_executable(${PROJECT_NAME} ${SOURCES_INC} ${SOURCES_CPP})


target_include_directories(${PROJECT_NAME} PRIVATE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  ${GTEST_INCLUDE_DIRS}
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  )

target_link_libraries(${PROJECT_NAME} PRIVATE
  phare_amr
  ${GTEST_LIBS})


add_phare_test(${PROJECT_NAME} ${CMAKE_CURRENT_BINARY_DIR})


