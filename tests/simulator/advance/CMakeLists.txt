cmake_minimum_required (VERSION 3.9)

project(test-simulator-advance)

add_subdirectory(_1d)

if(NOT lowResourceTests)

  add_subdirectory(_2d)

endif()

