
if (test AND ${PHARE_EXEC_LEVEL_MIN} GREATER 0) # 0 = no tests

  # each directory adds its own subdirectories
  add_subdirectory(tests)
  add_subdirectory(pyphare/pyphare_tests)

endif()
