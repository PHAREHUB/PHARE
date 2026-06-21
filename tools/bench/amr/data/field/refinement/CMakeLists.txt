cmake_minimum_required (VERSION 3.20.1)

project(phare_bench_field_refine)

add_phare_cpp_benchmark(10 ${PROJECT_NAME} bench_field_refinement ${CMAKE_CURRENT_BINARY_DIR})

