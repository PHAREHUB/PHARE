cmake_minimum_required (VERSION 3.9)

project(phare_python3)

pybind11_add_module(cpp cpp_simulator.cpp)
target_link_libraries(cpp PUBLIC phare_simulator)
target_compile_options(cpp PRIVATE ${PHARE_FLAGS} -DPHARE_HAS_HIGHFIVE=${PHARE_HAS_HIGHFIVE}) # pybind fails with Werror
set_target_properties(cpp
    PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/pybindlibs"
)
# this is on by default "pybind11_add_module" but can interfere with coverage so we disable it if coverage is enabled
set_property(TARGET cpp PROPERTY INTERPROCEDURAL_OPTIMIZATION ${PHARE_INTERPROCEDURAL_OPTIMIZATION})
