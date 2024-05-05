cmake_minimum_required (VERSION 3.20.1)

project(phare_bench_ion_updater)

add_phare_cpp_benchmark(11 ${PROJECT_NAME} bench_ion_updater ${CMAKE_CURRENT_BINARY_DIR})
