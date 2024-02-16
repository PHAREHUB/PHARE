

set(CPPDICT_SRCDIR ${CMAKE_CURRENT_SOURCE_DIR}/subprojects/cppdict)

if (NOT EXISTS ${CPPDICT_SRCDIR})
  execute_process(
    COMMAND ${Git} clone https://github.com/LaboratoryOfPlasmaPhysics/cppdict ${CPPDICT_SRCDIR} -b master --recursive --depth 10
  )
else()
  if(devMode)
    message("downloading latest cppdict updates")
    execute_process(COMMAND ${Git} pull origin master WORKING_DIRECTORY ${CPPDICT_SRCDIR})
  endif(devMode)
endif()
