cmake_minimum_required (VERSION 3.20.1)

project(phare_bench_field)

add_phare_cpp_benchmark(10 ${PROJECT_NAME} field_box ${CMAKE_CURRENT_BINARY_DIR})
