cmake_minimum_required (VERSION 3.9)

project(phare_real_bench)

pybind11_add_module(cpp_sim_2_1_4 sim/sim_2_1_4.cpp)
target_link_libraries(cpp_sim_2_1_4 PUBLIC phare_simulator)
target_compile_options(cpp_sim_2_1_4 PRIVATE ${PHARE_FLAGS} -DPHARE_HAS_HIGHFIVE=${PHARE_HAS_HIGHFIVE}) # pybind fails with Werror
set_target_properties(cpp_sim_2_1_4
    PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/pybindlibs"
)
# this is on by default "pybind11_add_module" but can interfere with coverage so we disable it if coverage is enabled
set_property(TARGET cpp_sim_2_1_4 PROPERTY INTERPROCEDURAL_OPTIMIZATION ${PHARE_INTERPROCEDURAL_OPTIMIZATION})

