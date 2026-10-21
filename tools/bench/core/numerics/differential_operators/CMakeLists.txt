cmake_minimum_required (VERSION 3.20.1)

project(phare_bench_differential_operators)

add_phare_cpp_benchmark(10 ${PROJECT_NAME} bench_differential_operators ${CMAKE_CURRENT_BINARY_DIR})
