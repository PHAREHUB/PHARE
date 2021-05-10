cmake_minimum_required (VERSION 3.9)

project(test-simulator-init)

add_subdirectory(_1d)

if(NOT lowResourceTests)

  add_subdirectory(_2d)

endif()

