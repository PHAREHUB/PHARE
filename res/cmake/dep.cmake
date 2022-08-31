
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
find_program(Git git)

# SAMRAI build option
include("${PHARE_PROJECT_DIR}/res/cmake/dep/samrai.cmake")
include("${PHARE_PROJECT_DIR}/res/cmake/dep/raja_umpire.cmake")


# caliper build option
#  enabled with -DCALIPER_ROOT=/path/to/caliper
#    or -DwithCaliper, which dowloads to subprojects dir
include("${PHARE_PROJECT_DIR}/res/cmake/dep/caliper.cmake")

# pybind
include("${PHARE_PROJECT_DIR}/res/cmake/dep/pybind.cmake")


# HighFive
include("${PHARE_PROJECT_DIR}/res/cmake/dep/highfive.cmake")

include("${PHARE_PROJECT_DIR}/res/cmake/dep/gtest.cmake")
include("${PHARE_PROJECT_DIR}/res/cmake/dep/gbench.cmake")
