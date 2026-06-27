cmake_minimum_required (VERSION 3.20.1)

project(phare_bench_boxes)

add_phare_cpp_benchmark(10 ${PROJECT_NAME} box_iterator ${CMAKE_CURRENT_BINARY_DIR})
add_phare_cpp_benchmark(10 ${PROJECT_NAME} box_span ${CMAKE_CURRENT_BINARY_DIR})
add_phare_cpp_benchmark(10 ${PROJECT_NAME} box_span_N_loop ${CMAKE_CURRENT_BINARY_DIR})

