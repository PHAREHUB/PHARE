cmake_minimum_required (VERSION 3.20.1)

project(phare_bench_utilities)

add_phare_cpp_benchmark(10 ${PROJECT_NAME} timestamps ${CMAKE_CURRENT_BINARY_DIR})
