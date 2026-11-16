if(NOT DEFINED SKBUILD OR NOT CMAKE_GENERATOR MATCHES "Ninja")
  return()
endif()

set(_phare_ninja_launcher "${CMAKE_BINARY_DIR}/phare_ninja")
if(NOT CMAKE_MAKE_PROGRAM STREQUAL _phare_ninja_launcher)
  set(PHARE_NINJA "${CMAKE_MAKE_PROGRAM}" CACHE FILEPATH "Ninja run behind the pip/uv progress line" FORCE)
endif()

file(CONFIGURE OUTPUT "${_phare_ninja_launcher}" CONTENT
  "#!/bin/sh\nexec \"@Python_EXECUTABLE@\" \"@PHARE_PROJECT_DIR@/res/cmake/ninja_progress.py\" \"@CMAKE_BINARY_DIR@\" \"@PHARE_NINJA@\" \"$@\"\n"
  @ONLY)
file(CHMOD "${_phare_ninja_launcher}"
  PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
set(CMAKE_MAKE_PROGRAM "${_phare_ninja_launcher}" CACHE FILEPATH "Program used to build from build.ninja files." FORCE)
