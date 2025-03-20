

function(phare_git_get_or_update name dir url branch)
  if (NOT EXISTS ${dir})
    message("cloning ${url} ${branch}" )
    execute_process(
      COMMAND ${Git} clone ${url} ${dir} -b ${branch} --recursive --depth 1 --shallow-submodules
    )
  else()
    if(devMode)
      message("downloading latest ${name} updates")
      execute_process(COMMAND ${Git} pull origin ${branch} WORKING_DIRECTORY ${dir})
    endif(devMode)
  endif()
endfunction(phare_git_get_or_update)

function(phare_github_get_or_update name dir repo branch)
  phare_git_get_or_update(${name} ${dir} https://github.com/${repo} ${branch})
endfunction(phare_github_get_or_update)

include("${PHARE_PROJECT_DIR}/res/cmake/dep/cppdict.cmake")
add_subdirectory(subprojects/cppdict)

# SAMRAI build option
include("${PHARE_PROJECT_DIR}/res/cmake/dep/samrai.cmake")


# caliper build option
#  enabled with -DCALIPER_ROOT=/path/to/caliper
#    or -DwithCaliper, which downloads to subprojects dir
include("${PHARE_PROJECT_DIR}/res/cmake/dep/caliper.cmake")

# pybind
include("${PHARE_PROJECT_DIR}/res/cmake/dep/pybind.cmake")


# HighFive
include("${PHARE_PROJECT_DIR}/res/cmake/dep/highfive.cmake")

# Phlop
include("${PHARE_PROJECT_DIR}/res/cmake/dep/phlop.cmake")

