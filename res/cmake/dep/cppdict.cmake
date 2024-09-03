

set(CPPDICT_SRCDIR ${CMAKE_CURRENT_SOURCE_DIR}/subprojects/cppdict)
set(CPPDICT_VERSION master)
phare_github_get_or_update(cppdict ${CPPDICT_SRCDIR} LaboratoryOfPlasmaPhysics/cppdict ${CPPDICT_VERSION})

