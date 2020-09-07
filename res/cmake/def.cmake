

set (PHARE_FLAGS ${PHARE_FLAGS})
set (PHARE_WERROR_FLAGS ${PHARE_FLAGS} ${PHARE_WERROR_FLAGS})
set (PHARE_PYTHONPATH "${CMAKE_BINARY_DIR}:${CMAKE_SOURCE_DIR}/pyphare")

# Link Time Optimisation flags - is disabled if coverage is enabled
set (PHARE_INTERPROCEDURAL_OPTIMIZATION FALSE)
if(withIPO)
  include(CheckIPOSupported)
  check_ipo_supported(RESULT PHARE_INTERPROCEDURAL_OPTIMIZATION OUTPUT error)
endif(withIPO)

if(devMode) # -DdevMode=ON
  # Having quotes on strings here has lead to quotes being added to the compile string, so avoid.

  set (_Werr ${PHARE_WERROR_FLAGS} -Wall -Wextra -pedantic -Werror -Wno-unused-variable -Wno-unused-parameter)

  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set (_Werr ${_Werr} -Wno-gnu-zero-variadic-macro-arguments)

  else() # !Clang
    set (_Werr ${_Werr} -Wno-class-memaccess -Wno-unused-but-set-variable -Wno-unused-but-set-parameter)

  endif() # clang

  set (PHARE_WERROR_FLAGS ${_Werr})

endif(devMode)

function(phare_sanitize_ san cflags )
  set(CMAKE_REQUIRED_FLAGS ${san})
  check_cxx_compiler_flag( ${san} ADDRESS_SANITIZER)
  if (${ADDRESS_SANITIZER})
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${san} ${cflags}" PARENT_SCOPE)
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${san}"  PARENT_SCOPE)
    set (CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${san}"  PARENT_SCOPE)
  else()
    message(FATAL_ERROR "Your compiler: ${CMAKE_CXX_COMPILER_ID} seems to not support ${san}")
  endif()
  unset(CMAKE_REQUIRED_FLAGS)
endfunction(phare_sanitize_)

if (asan)   # -Dasan=ON
  phare_sanitize_("-fsanitize=address" "-fno-omit-frame-pointer" )
endif(asan)

if (ubsan)  # -Dubsan=ON
  phare_sanitize_("-fsanitize=undefined" "" )
endif(ubsan)

# msan is not supported - it's not practical to configure - use valgrind
