cmake_minimum_required (VERSION 3.9)

project(test-simulator-init)

add_subdirectory(1d)

if(NOT lowResourceTests)

  add_subdirectory(2d)

endif()

