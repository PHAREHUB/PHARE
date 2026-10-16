cmake_minimum_required (VERSION 3.20.1)

project(phare_python3)

file(WRITE "${CMAKE_BINARY_DIR}/pybindlibs/__init__.py" "")

if(NOT PHARE_PERMUTATIONS)
  set(PHARE_PERMUTATIONS "${PHARE_PROJECT_DIR}/res/sim/all.txt")
endif()

if(NOT IS_ABSOLUTE "${PHARE_PERMUTATIONS}")
  set(PHARE_PERMUTATIONS "${PHARE_PROJECT_DIR}/${PHARE_PERMUTATIONS}")
endif()

if(NOT EXISTS "${PHARE_PERMUTATIONS}")
  message(FATAL_ERROR "PHARE_PERMUTATIONS file not found: ${PHARE_PERMUTATIONS}")
endif()

file(STRINGS ${PHARE_PERMUTATIONS} THIS_FILE)

if(DEFINED SKBUILD)
  unset(PHARE_PERMUTATIONS CACHE)
endif()

while(THIS_FILE)
    list(POP_FRONT THIS_FILE LINE)

    # Skip empty lines and comments
    string(STRIP "${LINE}" LINE)
    if("${LINE}" STREQUAL "" OR "${LINE}" MATCHES "^#" OR NOT "${LINE}" MATCHES "^[0-9]")
        continue()
    endif()

    string(REPLACE "," ";" LINE_LIST ${LINE})
    list(LENGTH LINE_LIST LINE_LENGTH)

    if(LINE_LENGTH EQUAL 5)
        # MHD-only shape: dim,recon,slope,riemann,hall
        # normalize to the 7-field SimOpts aggregate order: dim,0,0,...
        list(INSERT LINE_LIST 1 "0" "0")
        list(LENGTH LINE_LIST LINE_LENGTH)
    elseif(NOT LINE_LENGTH EQUAL 3 AND NOT LINE_LENGTH EQUAL 7)
        message(FATAL_ERROR "PHARE permutation line has unsupported arity (${LINE_LENGTH}): ${LINE}")
    endif()

    SET(_OFFSET 3)
    SET(_PREFIXES "ReconstructionType;SlopeLimiterType;RiemannSolverType")
    foreach(PREFIX IN LISTS _PREFIXES)
      if(LINE_LENGTH GREATER ${_OFFSET})
          list(TRANSFORM LINE_LIST PREPEND "${PREFIX}::" AT ${_OFFSET})
      endif()
      MATH(EXPR _OFFSET "${_OFFSET}+1")
    endforeach()

    string(REPLACE "," "_" OUT ${LINE})
    SET(MOD_STR "cpp_${OUT}")
    list(JOIN LINE_LIST "," LINE) # semicolon separated list to csv

    pybind11_add_module(${MOD_STR} cpp_simulator.cpp)
    target_link_libraries(${MOD_STR} PUBLIC phare_simulator)
    target_compile_options(${MOD_STR} PRIVATE
        ${PHARE_WERROR_FLAGS}
        -DPHARE_HAS_HIGHFIVE=${PHARE_HAS_HIGHFIVE}
        -DPHARE_SIM_ID=${OUT}
        "-DPHARE_SIM_STR=${LINE}"
    )
    set_target_properties(${MOD_STR}
        PROPERTIES
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/pybindlibs"
    )
    # this is on by default "pybind11_add_module" but can interfere with coverage so we disable it if coverage is enabled
    set_property(TARGET ${MOD_STR} PROPERTY INTERPROCEDURAL_OPTIMIZATION ${PHARE_INTERPROCEDURAL_OPTIMIZATION})
    set_property(TARGET ${MOD_STR} APPEND_STRING PROPERTY LINK_FLAGS " ${PHARE_LINK_FLAGS}")
    set_property(GLOBAL APPEND PROPERTY PHARE_PYTHON_MODULES ${MOD_STR})
    phare_limit_module_compile(${MOD_STR})
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
set_property(GLOBAL APPEND PROPERTY PHARE_PYTHON_MODULES cpp_etc)
phare_limit_module_compile(cpp_etc)
