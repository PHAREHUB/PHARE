
include("${PHARE_PROJECT_DIR}/res/cmake/dep/cppdict.cmake")
add_subdirectory(subprojects/cppdict)

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

# Phlop
include("${PHARE_PROJECT_DIR}/res/cmake/dep/phlop.cmake")

