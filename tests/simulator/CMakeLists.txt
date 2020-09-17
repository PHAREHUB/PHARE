cmake_minimum_required (VERSION 3.9)

project(test-simulator)

if(NOT ${PHARE_PROJECT_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  file(GLOB PYFILES "*.py")
  file(COPY ${PYFILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif()

if(HighFive)

  ## These test use dump diagnostics so require HighFive!
  add_python3_test(diagnostics     diagnostics.py         ${CMAKE_CURRENT_BINARY_DIR})
  add_python3_test(init-particles  test_initialization.py ${CMAKE_CURRENT_BINARY_DIR})
  add_python3_test(sim-refineboxes refinement_boxes.py    ${CMAKE_CURRENT_BINARY_DIR})

endif()

add_python3_test(data-wrangler        data_wrangler.py        ${CMAKE_CURRENT_BINARY_DIR})
add_python3_test(sim-refineParticlNbr refined_particle_nbr.py ${CMAKE_CURRENT_BINARY_DIR})

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.py ${CMAKE_CURRENT_BINARY_DIR}/config.py @ONLY)
