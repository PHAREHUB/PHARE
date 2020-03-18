cmake_minimum_required (VERSION 3.3)

project(test-simulator)

if(HighFive)

  if(NOT ${PHARE_PROJECT_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    file(COPY ${PROJECT_SOURCE_DIR}/py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  endif()

  add_python3_test(sim-refineboxes py/refinement_boxes.py ${CMAKE_CURRENT_BINARY_DIR})
  add_python3_test(data-wrangler   py/data_wrangler.py ${CMAKE_CURRENT_BINARY_DIR})

endif()
