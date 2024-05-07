
# public test functions
#
#   add_phare_cpp_benchmark(exec_level project file directory)
#    create $target ${project}_${file}
#    execute $target in $directory, with mpirun when -DtestMPI=ON if exec_level is active
#

if (bench)

  if (NOT DEFINED PHARE_BENCH_PROFILING)
    set(PHARE_BENCH_PROFILE 0)
  endif()

  if(DEFINED GBENCH_ROOT)
    set(GBENCH_ROOT ${GBENCH_ROOT} CACHE PATH "Path to googlebenchmark")
    find_package(benchmark REQUIRED)
    set(GBENCH_LIBS benchmark::benchmark)
  else()
    set(GBENCH_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/subprojects/googlebench)

    if (NOT EXISTS ${GBENCH_ROOT})
      execute_process(COMMAND ${Git} clone https://github.com/google/benchmark ${GBENCH_ROOT} --depth 1)
    endif()

    option(BENCHMARK_ENABLE_TESTING "Enable building the unit tests which depend on gtest" OFF)
    option(BENCHMARK_ENABLE_GTEST_TESTS "Enable building the unit tests which depend on gtest" OFF)
    add_subdirectory(subprojects/googlebench)
    set(GBENCH_LIBS "benchmark")

  endif()

  function(add_phare_cpp_benchmark_ exec_level target file directory)
    add_executable(${target} ${file})
    target_compile_options(${target} PRIVATE ${PHARE_WERROR_FLAGS} -DPHARE_HAS_HIGHFIVE=${PHARE_HAS_HIGHFIVE})
    set_property(TARGET ${target} PROPERTY INTERPROCEDURAL_OPTIMIZATION ${PHARE_INTERPROCEDURAL_OPTIMIZATION})
    target_include_directories(${target} PUBLIC ${PHARE_PROJECT_DIR}/subprojects/googlebench/include ${PHARE_PROJECT_DIR})
    add_phare_test(${target} ${directory}) # using this function means benchmarks can be run with MPI, not sure this is good.
  endfunction(add_phare_cpp_benchmark_)

  function(add_phare_cpp_benchmark exec_level project file directory)
    if(${exec_level} GREATER_EQUAL ${PHARE_EXEC_LEVEL_MIN} AND ${exec_level} LESS_EQUAL ${PHARE_EXEC_LEVEL_MAX})
      add_phare_cpp_benchmark_(${exec_level} ${project}_${file} ${file}.cpp ${directory})
      target_link_libraries(${project}_${file} PUBLIC ${GBENCH_LIBS} phare_simulator)
    endif()
  endfunction(add_phare_cpp_benchmark)


  add_subdirectory(tools/bench/core/data/particles)
  add_subdirectory(tools/bench/core/numerics/pusher)

  add_subdirectory(tools/bench/amr/data/particles)

  add_subdirectory(tools/bench/hi5)
  add_subdirectory(tools/bench/real)

endif()
