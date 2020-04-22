cmake_minimum_required(VERSION 3.0)
project(phare_initializer)




set(SOURCE_INC data_provider.h
               python_data_provider.h
               restart_data_provider.h)


set(SOURCE_CPP data_provider.cpp)

set(PYBIND11_CPP_STANDARD -std=c++17)

add_library(${PROJECT_NAME} SHARED ${SOURCE_CPP} ${SOURCE_INC})
target_compile_options(${PROJECT_NAME} PRIVATE ${PHARE_WERROR_FLAGS})


target_include_directories(${PROJECT_NAME}  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../subprojects>)


# this is the python module that wraps the add() method
pybind11_add_module(dictator dictator.cpp)
target_link_libraries(dictator PUBLIC phare_initializer pybind11::embed)
set_target_properties(dictator
    PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/pybindlibs"
)

#configure_file(${CMAKE_CURRENT_SOURCE_DIR}/__init__.py.in ${CMAKE_BINARY_DIR}/pyphare/pyphare/__init__.py @ONLY)
