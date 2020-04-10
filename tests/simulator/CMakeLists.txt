cmake_minimum_required (VERSION 3.3)

project(test-simulator)

if(HighFive)

  if(NOT ${PHARE_PROJECT_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    file(COPY ${PROJECT_SOURCE_DIR}/py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  endif()

  pybind11_add_module(test_simulator test_simulator_pybind.cpp)

  target_link_libraries(test_simulator PUBLIC
    phare_simulator
    ${GTEST_LIBS}
    ${HDF5_LIBRARIES})

  set_target_properties(test_simulator
      PROPERTIES
      LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
  )

  add_python3_test(sim-refineboxes  py/refinement_boxes.py ${CMAKE_CURRENT_BINARY_DIR})
  add_python3_test(sim-overlap-init py/overlap_init.py     ${CMAKE_CURRENT_BINARY_DIR})
  add_python3_test(sim-fluid-init   py/fluid_init.py       ${CMAKE_CURRENT_BINARY_DIR})
  add_python3_test(sim-em-init      py/em_init.py          ${CMAKE_CURRENT_BINARY_DIR})
  add_python3_test(data-wrangler    py/data_wrangler.py    ${CMAKE_CURRENT_BINARY_DIR})

  if(devMode)
    add_python3_test(sim-overlaps     py/overlaps.py         ${CMAKE_CURRENT_BINARY_DIR})
    if (NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/py/overlaps)
      add_custom_command(TARGET test_simulator
                         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/py
                         POST_BUILD
                         COMMAND git clone https://github.com/PHARCHIVE/hdf5_patch_overlap --depth 1 overlaps)
      add_custom_command(TARGET test_simulator
                         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/py/overlaps
                         POST_BUILD
                         COMMAND find . -not -path '*/\\.*' -type f -name "*.xz" | xargs unxz)
    endif() # clone/uncompress

  endif() # devMode

endif() # HighFive
