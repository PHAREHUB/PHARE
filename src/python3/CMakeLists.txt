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

if(devMode)

target_include_directories(cpp PRIVATE ${CMAKE_BINARY_DIR} ${CMAKE_CURRENT_BINARY_DIR})
target_compile_definitions(cpp PRIVATE PHARE_DEV_MODE=${PHARE_DEV_MODE})
add_custom_target(cpp.timestamp ALL COMMAND date +%s > cpp.timestamp WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
add_dependencies(cpp cpp.timestamp)

target_include_directories(data_wrangler PRIVATE ${CMAKE_BINARY_DIR} ${CMAKE_CURRENT_BINARY_DIR})
target_compile_definitions(data_wrangler PRIVATE PHARE_DEV_MODE=${PHARE_DEV_MODE})
add_custom_target(data_wrangler.timestamp ALL COMMAND date +%s > data_wrangler.timestamp WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
add_dependencies(data_wrangler data_wrangler.timestamp)

endif() # devMode
