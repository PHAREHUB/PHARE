
if (bench)

  if(DEFINED GBENCH_ROOT)
    set(GBENCH_ROOT ${GBENCH_ROOT} CACHE PATH "Path to googlebenchmark")
    find_package(benchmark REQUIRED)
    set(GBENCH_LIBS benchmark::benchmark)
  else()
    set(GBENCH_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/subprojects/googlebench)

    if (NOT EXISTS ${GBENCH_ROOT})
      execute_process(COMMAND ${Git} clone https://github.com/google/benchmark ${GBENCH_ROOT} --depth 1)
    else()
      if(devMode)
        message("downloading latest googlebench updates")
        execute_process(COMMAND ${Git} pull origin master WORKING_DIRECTORY ${GBENCH_ROOT})
      endif(devMode)
    endif()

    option(BENCHMARK_ENABLE_TESTING "Enable building the unit tests which depend on gtest" OFF)
    option(BENCHMARK_ENABLE_GTEST_TESTS "Enable building the unit tests which depend on gtest" OFF)
    add_subdirectory(subprojects/googlebench)
    set(GBENCH_LIBS "benchmark")

  endif()

endif()
