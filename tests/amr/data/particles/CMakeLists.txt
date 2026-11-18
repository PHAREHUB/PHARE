add_subdirectory(copy)
add_subdirectory(copy_overlap)
add_subdirectory(stream_pack)
add_subdirectory(refine)

if(HighFive) # the test needs SAMRAI with HDF5, expected with HighFive
  add_subdirectory(restart)
endif()
