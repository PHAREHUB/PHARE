

if(NOT DEFINED PHARE_CPPDICT_VERSION)
  SET(PHARE_CPPDICT_VERSION "master")
endif()

set(CPPDICT_SRCDIR ${CMAKE_CURRENT_SOURCE_DIR}/subprojects/cppdict)
phare_github_get_or_update(cppdict ${CPPDICT_SRCDIR} LaboratoryOfPlasmaPhysics/cppdict ${PHARE_CPPDICT_VERSION})
