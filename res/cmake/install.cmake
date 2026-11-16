if(APPLE)
  set(_phare_origin "@loader_path")
else()
  set(_phare_origin "$ORIGIN")
endif()

get_property(_phare_python_modules GLOBAL PROPERTY PHARE_PYTHON_MODULES)
foreach(_mod IN LISTS _phare_python_modules)
  set_target_properties(${_mod} PROPERTIES
    INSTALL_RPATH "${_phare_origin}/lib"
    INSTALL_RPATH_USE_LINK_PATH TRUE
  )
  install(TARGETS ${_mod} LIBRARY DESTINATION pybindlibs COMPONENT pyphare)
endforeach()
install(FILES "${CMAKE_BINARY_DIR}/pybindlibs/__init__.py" DESTINATION pybindlibs COMPONENT pyphare)

set(_phare_candidate_libs
  phare_core phare_mpi phare_amr phare_initializer
  SAMRAI_algs SAMRAI_appu SAMRAI_geom SAMRAI_hier SAMRAI_math
  SAMRAI_mesh SAMRAI_pdat SAMRAI_solv SAMRAI_tbox SAMRAI_xfer
)
foreach(_lib IN LISTS _phare_candidate_libs)
  if(NOT TARGET ${_lib})
    continue()
  endif()
  get_target_property(_type ${_lib} TYPE)
  get_target_property(_imported ${_lib} IMPORTED)
  if(NOT _type STREQUAL "SHARED_LIBRARY" OR _imported)
    continue()
  endif()
  set_target_properties(${_lib} PROPERTIES
    INSTALL_RPATH "${_phare_origin}"
    INSTALL_RPATH_USE_LINK_PATH TRUE
  )
  install(TARGETS ${_lib} LIBRARY DESTINATION pybindlibs/lib COMPONENT pyphare)
endforeach()
