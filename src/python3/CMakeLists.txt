cmake_minimum_required (VERSION 3.20.1)

project(phare_python3)

pybind11_add_module(pyMHD pyMHD.cpp)
target_include_directories(pyMHD 
    PRIVATE 
    ${CMAKE_SOURCE_DIR}
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../subprojects>
)
target_link_libraries(pyMHD PUBLIC phare_simulator)
target_compile_options(pyMHD PRIVATE ${PHARE_WERROR_FLAGS} -DPHARE_HAS_HIGHFIVE=${PHARE_HAS_HIGHFIVE})
set_target_properties(pyMHD
    PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/pybindlibs"
)

pybind11_add_module(cpp cpp_simulator.cpp)
target_link_libraries(cpp PUBLIC phare_simulator)
target_compile_options(cpp PRIVATE ${PHARE_WERROR_FLAGS} -DPHARE_HAS_HIGHFIVE=${PHARE_HAS_HIGHFIVE})
set_target_properties(cpp
    PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/pybindlibs"
)
# this is on by default "pybind11_add_module" but can interfere with coverage so we disable it if coverage is enabled
set_property(TARGET cpp PROPERTY INTERPROCEDURAL_OPTIMIZATION ${PHARE_INTERPROCEDURAL_OPTIMIZATION})
set_property(TARGET cpp APPEND_STRING PROPERTY LINK_FLAGS " ${PHARE_LINK_FLAGS}")

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
  pybind11_add_module(cpp_dbg cpp_simulator.cpp)
  target_link_libraries(cpp_dbg PUBLIC phare_simulator)
  target_compile_options(cpp_dbg PRIVATE ${PHARE_WERROR_FLAGS} -DPHARE_HAS_HIGHFIVE=${PHARE_HAS_HIGHFIVE} -DPHARE_DIAG_DOUBLES=1 -DPHARE_CPP_MOD_NAME=cpp_dbg)
  set_target_properties(cpp_dbg
      PROPERTIES
      LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/pybindlibs"
  )
  set_property(TARGET cpp_dbg PROPERTY INTERPROCEDURAL_OPTIMIZATION ${PHARE_INTERPROCEDURAL_OPTIMIZATION})
  set_property(TARGET cpp_dbg APPEND_STRING PROPERTY LINK_FLAGS " ${PHARE_LINK_FLAGS}")
endif (CMAKE_BUILD_TYPE STREQUAL "Debug")


pybind11_add_module(cpp_etc cpp_etc.cpp)
target_compile_options(cpp_etc PRIVATE ${PHARE_WERROR_FLAGS} -DPHARE_HAS_HIGHFIVE=${PHARE_HAS_HIGHFIVE})
target_link_libraries(cpp_etc PUBLIC phare_amr phare_diagnostic)
set_target_properties(cpp_etc
    PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/pybindlibs"
)
# this is on by default "pybind11_add_module" but can interfere with coverage so we disable it if coverage is enabled
set_property(TARGET cpp_etc PROPERTY INTERPROCEDURAL_OPTIMIZATION ${PHARE_INTERPROCEDURAL_OPTIMIZATION})
set_property(TARGET cpp_etc APPEND_STRING PROPERTY LINK_FLAGS " ${PHARE_LINK_FLAGS}")
