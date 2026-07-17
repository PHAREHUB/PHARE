


set(MAGIC_ENUM_DIR ${CMAKE_CURRENT_SOURCE_DIR}/subprojects/menum)
if (NOT EXISTS ${MAGIC_ENUM_DIR})
   # menum&subprojects/menum(https://github.com/Neargye/magic_enum)
  execute_process(COMMAND ${Git} clone https://github.com/Neargye/magic_enum ${MAGIC_ENUM_DIR})
else()
  if(devMode)
    message("downloading latest magic enum updates")
    execute_process(COMMAND ${Git} pull origin master WORKING_DIRECTORY ${MAGIC_ENUM_DIR})
  endif(devMode)
endif()

add_subdirectory(subprojects/menum)
