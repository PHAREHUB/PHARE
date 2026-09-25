
if (test AND coverage)

  # LTO disabled for coverage builds
  set (PHARE_INTERPROCEDURAL_OPTIMIZATION FALSE)

  set (_Fvr " -fprofile-arcs -ftest-coverage")

  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg -DHAVE_EXECINFO_H -g3 -O0 ${_Fvr}")
  set(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS}  ${_Fvr}")

  add_custom_target(build-time-make-directory ALL
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/coverage)

  set (_Gcvr gcovr --exclude=.*subprojects.* --exclude=.*tests.* --exclude=/usr/include/.* )
  set (_Gcvr ${_Gcvr} --object-directory ${CMAKE_BINARY_DIR} -r ${CMAKE_SOURCE_DIR})

  # hot functions (e.g. Field::operator()) legitimately exceed gcovr's default suspicious
  #  hits threshold (2^32) over the test suite. Raise it rather than disable it, so the
  #  near 2^64 counts of https://gcc.gnu.org/bugzilla/show_bug.cgi?id=68080 are still caught.
  #  The option only exists in newer gcovr versions.
  execute_process(COMMAND gcovr --help OUTPUT_VARIABLE _Gcvr_help ERROR_QUIET)
  if (_Gcvr_help MATCHES "--gcov-suspicious-hits-threshold")
    set (_Gcvr ${_Gcvr} --gcov-suspicious-hits-threshold 281474976710656) # 2^48
  endif()

  # one pass over the coverage data for both reports, always regenerated
  add_custom_target(gcovr
    COMMAND ${_Gcvr}
            --html-details ${CMAKE_CURRENT_BINARY_DIR}/coverage/index.html
            --xml ${CMAKE_CURRENT_BINARY_DIR}/coverage/coverage.xml
    DEPENDS build-time-make-directory
    VERBATIM
  )

  if(APPLE)
    set(OPPEN_CMD open)
  elseif(UNIX)
    set(OPPEN_CMD xdg-open)
  endif(APPLE)

  add_custom_target(show_coverage
    COMMAND ${OPPEN_CMD} ${CMAKE_CURRENT_BINARY_DIR}/coverage/index.html
    DEPENDS gcovr
  )

ENDIF(test AND coverage)
