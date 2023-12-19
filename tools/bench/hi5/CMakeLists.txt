cmake_minimum_required (VERSION 3.20.1)

project(phare_bench_h5)

if(HighFive)

  add_phare_cpp_benchmark(11 ${PROJECT_NAME} write_particles ${CMAKE_CURRENT_BINARY_DIR})

endif()
