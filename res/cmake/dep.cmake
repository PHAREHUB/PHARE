


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


set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
find_program(Git git REQUIRED)


# cppdict
include("${PHARE_PROJECT_DIR}/res/cmake/dep/cppdict.cmake")

# HighFive
include("${PHARE_PROJECT_DIR}/res/cmake/dep/highfive.cmake")

# SAMRAI/etc
include("${PHARE_PROJECT_DIR}/res/cmake/dep/raja_umpire.cmake")
include("${PHARE_PROJECT_DIR}/res/cmake/dep/samrai.cmake")


# caliper
#  enabled with -DCALIPER_ROOT=/path/to/caliper
#    or -DwithCaliper, which downloads to subprojects dir
include("${PHARE_PROJECT_DIR}/res/cmake/dep/caliper.cmake")

# pybind
include("${PHARE_PROJECT_DIR}/res/cmake/dep/pybind.cmake")


# Phlop - enabled with -DwithPhlop
include("${PHARE_PROJECT_DIR}/res/cmake/dep/phlop.cmake")

# google test
include("${PHARE_PROJECT_DIR}/res/cmake/dep/gtest.cmake")
# google benchmark
include("${PHARE_PROJECT_DIR}/res/cmake/dep/gbench.cmake")

include("${PHARE_PROJECT_DIR}/res/cmake/dep/magicenum.cmake")

# KokkosTools
include("${PHARE_PROJECT_DIR}/res/cmake/dep/kokkos_tools.cmake")

include("${PHARE_PROJECT_DIR}/res/cmake/dep/mkn.cmake")

phare_github_get_or_update(bstp ${PHARE_PROJECT_DIR}/subprojects/bstp bshoshany/thread-pool master)
include_directories(${PHARE_PROJECT_DIR}/subprojects/bstp/include)
