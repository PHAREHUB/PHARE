cmake_minimum_required (VERSION 3.3)

project(phare_python3)

pybind11_add_module(cpp cpp_simulator.cpp)
target_link_libraries(cpp PUBLIC phare_simulator)
set_target_properties(cpp
    PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/phare"
)

pybind11_add_module(data_wrangler data_wrangler.cpp)
target_link_libraries(data_wrangler PUBLIC phare_simulator)
set_target_properties(data_wrangler
    PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/phare/data"
)

