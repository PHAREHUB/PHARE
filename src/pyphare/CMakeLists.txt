
cmake_minimum_required(VERSION 3.0)
project(pyphare)

set(PYBIND11_CPP_STANDARD -std=c++17)

pybind11_add_module(pyphare pyphare.cpp)

target_include_directories(${PROJECT_NAME}  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include/phare/initializer>)

# this is the python module that wraps the add() method
target_link_libraries(pyphare PUBLIC initializer phare_amr phare_core phare_diagnostic pybind11::embed)
