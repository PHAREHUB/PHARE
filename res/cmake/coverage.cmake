
if (test)

  IF(coverage)
      set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -g -O0 -Wall -W -Wshadow -Wunused-variable -Wunused-parameter -Wunused-function -Wunused -Wno-system-headers -Wno-deprecated -Woverloaded-virtual -Wwrite-strings -fprofile-arcs -ftest-coverage")
      set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -g -O0 -Wall -W -Wshadow -Wunused-variable \
          -Wunused-parameter -Wunused-function -Wunused -Wno-system-headers \
          -Wno-deprecated -Woverloaded-virtual -Wwrite-strings -fprofile-arcs -ftest-coverage")
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fprofile-arcs -ftest-coverage")

      add_custom_target(build-time-make-directory ALL
          COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/coverage)

      add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/coverage/index.html
          COMMAND gcovr --exclude='.*subprojects.*' --exclude='.*tests.*' --exclude='/usr/include/.*' --object-directory ${CMAKE_BINARY_DIR}  -r ${CMAKE_SOURCE_DIR} --html  --html-details -o ${CMAKE_CURRENT_BINARY_DIR}/coverage/index.html
          )

    add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/coverage/coverage.xml
          COMMAND gcovr --exclude='.*subprojects.*' --exclude='.*tests.*' --exclude='/usr/include/.*' --object-directory ${CMAKE_BINARY_DIR}  -r ${CMAKE_SOURCE_DIR} --xml -o ${CMAKE_CURRENT_BINARY_DIR}/coverage/coverage.xml
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
  ENDIF(coverage)

ENDIF(test)
