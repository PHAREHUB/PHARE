
cmake_minimum_required (VERSION 3.20.1)
project(phare-exe)

add_executable(phare-exe ${SOURCES_INC} phare.cpp)
target_compile_options(phare-exe PRIVATE ${PHARE_WERROR_FLAGS} -DPHARE_HAS_HIGHFIVE=${PHARE_HAS_HIGHFIVE})
target_link_libraries(phare-exe PUBLIC
    ${PHARE_BASE_LIBS}
    phare_simulator
    pybind11::embed
)
set_property(TARGET ${PROJECT_NAME} APPEND_STRING PROPERTY LINK_FLAGS " ${PHARE_LINK_FLAGS}")

