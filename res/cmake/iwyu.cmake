if (withIWYU)

  # included after the dependencies, their targets don't get it, iwyu.py also skips anything
  # outside the --only directories (subprojects added later such as google benchmark)
  find_program(IWYU NAMES include-what-you-use iwyu REQUIRED)
  find_program(IWYU_FIX_INCLUDE NAMES fix_include fix_includes.py iwyu-fix-includes REQUIRED)
  # the version PHARE formats with
  find_program(IWYU_CLANG_FORMAT NAMES clang-format-20 clang-format REQUIRED)

  # the diff of what iwyu suggests, empty if nothing: ${CMAKE_BINARY_DIR}/iwyu.diff
  set(CMAKE_CXX_INCLUDE_WHAT_YOU_USE
    ${Python_EXECUTABLE} ${PHARE_PROJECT_DIR}/res/cmake/iwyu.py
      --iwyu ${IWYU}
      --fix_include ${IWYU_FIX_INCLUDE}
      # keeps the existing include order, what's added is then placed by .clang-format
      "--fix_args=--nocomments --nosafe_headers --noreorder"
      --clang_format ${IWYU_CLANG_FORMAT}
      --project ${PHARE_PROJECT_DIR}
      --build ${CMAKE_BINARY_DIR}
      --compile_commands ${CMAKE_BINARY_DIR}/compile_commands.json
      --only src
      --only tests
      # not tools/config: a CMake project of its own, run before this one
      --only tools/bench
      # headers on their own, with the flags of a translation unit that has every include path
      --headers_with src/python3/cpp_simulator.cpp
      --jobs 4
      --
      -Xiwyu --mapping_file=${PHARE_PROJECT_DIR}/res/iwyu.imp
      # includes, not forward declarations: templates need the definition where instantiated
      -Xiwyu --no_fwd_decls
  )

endif()
