

if(DEFINED WITH_UMPIRE AND DEFINED WITH_RAJA)
  find_package(RAJA PATHS ${WITH_RAJA} REQUIRED)
  if (NOT RAJA_FOUND)
    message(FATAL_ERROR "RAJA is needed")
  endif()

  find_package(umpire PATHS ${WITH_UMPIRE} REQUIRED)
  if (NOT umpire_FOUND)
    message(FATAL_ERROR "umpire is needed")
  endif()

  set (PHARE_BASE_LIBS ${PHARE_BASE_LIBS} RAJA umpire cudart)
  set (PHARE_FLAGS ${PHARE_FLAGS} -DHAVE_RAJA=1 -DHAVE_UMPIRE=1)

  set (THRUST_DIR ${PHARE_PROJECT_DIR}/subprojects/thrust)
  if (NOT EXISTS ${THRUST_DIR})
    execute_process(COMMAND ${Git} clone https://github.com/NVIDIA/thrust -b main ${THRUST_DIR} --depth 5 --shallow-submodules --recursive WORKING_DIRECTORY ${GTEST_ROOT})
  endif()
  include_directories(${THRUST_DIR})

endif()
