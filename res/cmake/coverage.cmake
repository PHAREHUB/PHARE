
if (test AND coverage)

  # LTO disabled for coverage builds
  set (PHARE_INTERPROCEDURAL_OPTIMIZATION FALSE)

  set (_Cvr " -g -O0 -Wall -W -Wshadow -Wunused-variable")
  set (_Cvr " ${_Cvr} -Wunused-parameter -Wunused-function -Wunused -Wno-system-headers")
  set (_Cvr " ${_Cvr} -Wno-deprecated -Woverloaded-virtual -Wwrite-strings")
  set (_Fvr " -fprofile-arcs -ftest-coverage")

  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg -DHAVE_EXECINFO_H -g3 ")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${_Cvr} ${_Fvr}")
  set(CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_C_FLAGS_DEBUG} ${_Cvr} ${_Fvr}")
  set(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS}  ${_Fvr}")

  add_custom_target(build-time-make-directory ALL
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/coverage)

  set (_Gcvr gcovr --exclude='.*subprojects.*' --exclude='.*tests.*' --exclude='/usr/include/.*' )
  set (_Gcvr ${_Gcvr} --object-directory ${CMAKE_BINARY_DIR} -r ${CMAKE_SOURCE_DIR})

  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/coverage/index.html
    COMMAND ${_Gcvr} --html --html-details -o ${CMAKE_CURRENT_BINARY_DIR}/coverage/index.html
  )

  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/coverage/coverage.xml
    COMMAND ${_Gcvr} --xml -o ${CMAKE_CURRENT_BINARY_DIR}/coverage/coverage.xml
  )

  add_custom_target(gcovr
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/coverage/index.html gcovr ${CMAKE_CURRENT_BINARY_DIR}/coverage/coverage.xml
  )

  if(APPLE)
    set(OPPEN_CMD open)
  elseif(UNIX)
    set(OPPEN_CMD xdg-open)
  endif(APPLE)

  add_custom_target(show_coverage
    COMMAND ${OPPEN_CMD} ${CMAKE_CURRENT_BINARY_DIR}/coverage/index.html
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/coverage/index.html gcovr
  )

ENDIF(test AND coverage)
