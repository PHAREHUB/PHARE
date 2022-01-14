cmake_minimum_required(VERSION 3.9)
project(phare_initializer)




set(SOURCE_INC data_provider.hpp
               python_data_provider.hpp
               restart_data_provider.hpp)


set(SOURCE_CPP data_provider.cpp)

set(PYBIND11_CPP_STANDARD -std=c++17)

add_library(${PROJECT_NAME} SHARED ${SOURCE_CPP} ${SOURCE_INC})
target_compile_options(${PROJECT_NAME} PRIVATE ${PHARE_WERROR_FLAGS})
set_property(TARGET ${PROJECT_NAME} PROPERTY INTERPROCEDURAL_OPTIMIZATION ${PHARE_INTERPROCEDURAL_OPTIMIZATION})


target_include_directories(${PROJECT_NAME}  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../subprojects>)


# this is the python module that wraps the add() method
pybind11_add_module(dictator dictator.cpp)
target_link_libraries(dictator PUBLIC phare_initializer pybind11::embed)
set_property(TARGET dictator PROPERTY INTERPROCEDURAL_OPTIMIZATION ${PHARE_INTERPROCEDURAL_OPTIMIZATION})
set_target_properties(dictator
    PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/pybindlibs"
)
