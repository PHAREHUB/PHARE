cmake_minimum_required (VERSION 3.9)

project(test-simulator-init)

if(NOT ${PHARE_PROJECT_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  file(GLOB PYFILES "*.py")
  file(COPY ${PYFILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif()

if(HighFive)
  ## These test use dump diagnostics so require HighFive!

  add_python3_test(init-1d-fields     test_fields_init_1d.py     ${CMAKE_CURRENT_BINARY_DIR})
  add_python3_test(init-1d-particles  test_particles_init_1d.py  ${CMAKE_CURRENT_BINARY_DIR})

  if(NOT lowResourceTests)

    add_python3_test(init-2d-fields     test_fields_init_2d.py     ${CMAKE_CURRENT_BINARY_DIR})
    add_python3_test(init-2d-particles  test_particles_init_2d.py  ${CMAKE_CURRENT_BINARY_DIR})

  endif()

endif()
