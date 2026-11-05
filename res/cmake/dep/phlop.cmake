
if (withPhlop)

  # pinned: the python package (pyproject.toml "timing") must read what these headers write
  if(NOT DEFINED PHARE_PHLOP_VERSION)
    SET(PHARE_PHLOP_VERSION "v0.0.45")
  endif()

  set(PHLOP_SRCDIR ${CMAKE_CURRENT_SOURCE_DIR}/subprojects/phlop)
  phare_github_get_or_update(phlop ${PHLOP_SRCDIR} PhilipDeegan/phlop ${PHARE_PHLOP_VERSION})
  include_directories(${PHLOP_SRCDIR}/inc)

endif(withPhlop)
