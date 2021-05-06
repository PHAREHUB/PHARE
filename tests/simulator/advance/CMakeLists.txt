cmake_minimum_required (VERSION 3.9)

project(test-simulator-advance)

add_subdirectory(1d)

if(NOT lowResourceTests)

  add_subdirectory(2d)

endif()

