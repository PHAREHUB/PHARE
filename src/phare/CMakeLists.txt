
cmake_minimum_required(VERSION 3.0)
project(phare-exe)

add_executable(phare-exe phare.cpp)

target_include_directories(initializer PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)


target_link_libraries(phare-exe PUBLIC phare_core phare_samrai_interface initializer)


# this is the script that PHARE runs to prepare pythonpath in case pyphare is
# ran in a virtualenv
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/setpythonpath.py.in ${CMAKE_CURRENT_BINARY_DIR}/setpythonpath.py @ONLY)
