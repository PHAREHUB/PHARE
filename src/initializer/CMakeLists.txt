cmake_minimum_required(VERSION 3.0)
project(phare_initializer)




set(SOURCE_INC data_provider.h
               python_data_provider.h
               restart_data_provider.h)


set(SOURCE_CPP python_data_provider.cpp data_provider.cpp)

set(PYBIND11_CPP_STANDARD -std=c++17)

add_library(${PROJECT_NAME} SHARED ${SOURCE_CPP} ${SOURCE_INC})
target_compile_options(${PROJECT_NAME} PRIVATE ${PHARE_WERROR_FLAGS})



target_include_directories(${PROJECT_NAME}  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../subprojects>)

target_link_libraries(${PROJECT_NAME} PUBLIC pybind11::embed)


# this is the python module that wraps the add() method
pybind11_add_module(pyphare pyphare.cpp)
target_link_libraries(pyphare PUBLIC phare_initializer)
set_target_properties(pyphare
    PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/phare"
)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/__init__.py.in ${CMAKE_BINARY_DIR}/phare/__init__.py @ONLY)
