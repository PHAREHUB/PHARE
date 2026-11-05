
if (test AND ${PHARE_EXEC_LEVEL_MIN} GREATER 0) # 0 = no tests

  # each directory adds its own subdirectories
  add_subdirectory(tests)
  add_subdirectory(pyphare/pyphare_tests)

  file(COPY "res/eg" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/res")
  add_subdirectory(res/eg)

endif()
