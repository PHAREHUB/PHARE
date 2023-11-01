

# SAMRAI build option
include("${PHARE_PROJECT_DIR}/res/cmake/dep/samrai.cmake")


# caliper build option
#  enabled with -DCALIPER_ROOT=/path/to/caliper
#    or -DwithCaliper, which dowloads to subprojects dir
include("${PHARE_PROJECT_DIR}/res/cmake/dep/caliper.cmake")

# pybind
include("${PHARE_PROJECT_DIR}/res/cmake/dep/pybind.cmake")


# HighFive
include("${PHARE_PROJECT_DIR}/res/cmake/dep/highfive.cmake")


add_subdirectory(subprojects/cppdict)
