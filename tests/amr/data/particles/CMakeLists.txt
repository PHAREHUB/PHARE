add_subdirectory(copy)
add_subdirectory(copy_overlap)
add_subdirectory(stream_pack)
add_subdirectory(refine)

function(_add_amr_particles_test src_name)

    add_executable(${src_name} ${src_name}.cpp)

    target_include_directories(${src_name} PRIVATE
      $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
      ${GTEST_INCLUDE_DIRS}
    )

    target_link_directories(${src_name} PRIVATE
      ${HDF5_LIBRARY_PATH}
    )

    target_link_libraries(${src_name} PRIVATE
      phare_amr
      ${GTEST_LIBS}
    )

    add_phare_test(${src_name} ${CMAKE_CURRENT_BINARY_DIR})

endfunction(_add_amr_particles_test)

file(COPY "test_particles_schedules_inputs" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

_add_amr_particles_test(test_particles_data)
_add_amr_particles_test(test_particles_schedules)

