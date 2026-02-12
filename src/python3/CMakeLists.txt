cmake_minimum_required (VERSION 3.20.1)

project(phare_python3)

if(NOT PHARE_PERMUTATIONS)
  set(PHARE_PERMUTATIONS "${PHARE_PROJECT_DIR}/res/sim/all.txt")
endif()

file(STRINGS ${PHARE_PERMUTATIONS} THIS_FILE)
while(THIS_FILE)
    list(POP_FRONT THIS_FILE LINE)
    string(REPLACE "," "_" OUT ${LINE})
    SET(MOD_STR "cpp_${OUT}")
    pybind11_add_module(${MOD_STR} cpp_simulator.cpp)
    target_link_libraries(${MOD_STR} PUBLIC phare_simulator )
    target_compile_options(${MOD_STR} PRIVATE
        ${PHARE_WERROR_FLAGS} -DPHARE_HAS_HIGHFIVE=${PHARE_HAS_HIGHFIVE} -DPHARE_SIM_ID=${OUT} -DPHARE_SIM_STR=${LINE}
    )
    set_target_properties(${MOD_STR}
        PROPERTIES
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/pybindlibs"
    )
    # this is on by default "pybind11_add_module" but can interfere with coverage so we disable it if coverage is enabled
    set_property(TARGET ${MOD_STR} PROPERTY INTERPROCEDURAL_OPTIMIZATION ${PHARE_INTERPROCEDURAL_OPTIMIZATION})
    set_property(TARGET ${MOD_STR} APPEND_STRING PROPERTY LINK_FLAGS " ${PHARE_LINK_FLAGS}")

endwhile()


# cpp_etc == everything not simulator related!
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
