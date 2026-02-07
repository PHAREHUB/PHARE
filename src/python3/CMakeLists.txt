cmake_minimum_required (VERSION 3.20.1)
project(phare_python3)

if(NOT PHARE_PERMUTATIONS)
  set(PHARE_PERMUTATIONS "${PHARE_PROJECT_DIR}/res/sim/all.txt")
endif()

file(STRINGS ${PHARE_PERMUTATIONS} THIS_FILE)

while(THIS_FILE)
    list(POP_FRONT THIS_FILE LINE)

    # Skip empty lines and comments
    string(STRIP "${LINE}" LINE)
    if("${LINE}" STREQUAL "" OR "${LINE}" MATCHES "^#")
        continue()
    endif()

    string(REPLACE "," ";" LINE_LIST ${LINE})
    list(LENGTH LINE_LIST LINE_LENGTH)

    list(GET LINE_LIST 0 DIM)
    list(GET LINE_LIST 1 INTERP)
    list(GET LINE_LIST 2 NBREF)

    if(LINE_LENGTH EQUAL 10)
        list(GET LINE_LIST 3 TIME_INT)
        list(GET LINE_LIST 4 RECON)
        list(GET LINE_LIST 5 SLOPE_LIM)
        list(GET LINE_LIST 6 RIEMANN)
        list(GET LINE_LIST 7 HALL)
        list(GET LINE_LIST 8 RESIST)
        list(GET LINE_LIST 9 HYPER_RESIST)

        set(PHARE_SIM_STR_EXPANDED
            "${DIM}, ${INTERP}, ${NBREF}, PHARE::MHDOpts::TimeIntegratorType::${TIME_INT}, PHARE::MHDOpts::ReconstructionType::${RECON}, PHARE::MHDOpts::SlopeLimiterType::${SLOPE_LIM}, PHARE::MHDOpts::RiemannSolverType::${RIEMANN}, ${HALL}, ${RESIST}, ${HYPER_RESIST}")

    # Hybrid Only
    elseif(LINE_LENGTH EQUAL 3)
        set(PHARE_SIM_STR_EXPANDED "${DIM}, ${INTERP}, ${NBREF}")
    else()
        message(FATAL_ERROR "Invalid line format in ${PHARE_PERMUTATIONS}: ${LINE}
Expected formats:
  - Hybrid Only: dim,interp,nbRef
  - With MHD: dim,interp,nbRef,TimeInt,Recon,SlopeLim,Riemann,Hall,Resist,HyperResist")
    endif()

    string(REPLACE "," "_" OUT ${LINE})
    SET(MOD_STR "cpp_${OUT}")

    pybind11_add_module(${MOD_STR} cpp_simulator.cpp)
    target_link_libraries(${MOD_STR} PUBLIC phare_simulator)
    target_compile_options(${MOD_STR} PRIVATE
        ${PHARE_WERROR_FLAGS}
        -DPHARE_HAS_HIGHFIVE=${PHARE_HAS_HIGHFIVE}
        -DPHARE_SIM_ID=${OUT}
        "-DPHARE_SIM_STR=${PHARE_SIM_STR_EXPANDED}"
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
