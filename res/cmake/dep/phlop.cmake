
if (withPhlop)

  set(PHLOP_SRCDIR ${CMAKE_CURRENT_SOURCE_DIR}/subprojects/phlop)
  phare_github_get_or_update(phlop ${PHLOP_SRCDIR} PhilipDeegan/phlop master)
  include_directories(${PHLOP_SRCDIR}/inc)

endif(withPhlop)
