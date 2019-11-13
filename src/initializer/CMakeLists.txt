
cmake_minimum_required(VERSION 3.0)
project(initializer)




set(SOURCE_INC data_provider.h
               python_data_provider.h
               restart_data_provider.h)


set(SOURCE_CPP python_data_provider.cpp data_provider.cpp)

set(PYBIND11_CPP_STANDARD -std=c++17)

add_library(initializer SHARED ${SOURCE_CPP} ${SOURCE_INC})
target_compile_options(initializer PRIVATE ${PHARE_WERROR_FLAGS})



target_include_directories(${PROJECT_NAME}  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../subprojects>
  $<INSTALL_INTERFACE:include/phare/initializer>)

target_link_libraries(initializer PUBLIC pybind11::embed)


# this is the python module that wraps the add() method
pybind11_add_module(pyphare pyphare.cpp)
target_link_libraries(pyphare PUBLIC initializer)


