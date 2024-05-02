
#  TO activate caliper with PHARE
#   1. configure cmake -DwithCaliper=ON
#   2. export CALI_CONFIG=runtime-report
#   or for inclusive timings per function (makes more sense IMO)
#   3. export CALI_CONFIG=runtime-report,calc.inclusive
##

if(DEFINED CALIPER_ROOT)
  set(withCaliper ON)
endif()

if (withCaliper)
  if(DEFINED CALIPER_ROOT)
    find_package(caliper PATHS ${CALIPER_ROOT} REQUIRED)
    include_directories("${caliper_INCLUDE_DIR}")
  else()
    set(CALIPER_SRCDIR ${CMAKE_CURRENT_SOURCE_DIR}/subprojects/caliper)
    set(CALIPER_BIN ${CMAKE_CURRENT_BINARY_DIR}/subprojects/caliper)

    if (NOT EXISTS ${CALIPER_SRCDIR})
      execute_process(
        COMMAND ${Git} clone https://github.com/LLNL/Caliper ${CALIPER_SRCDIR} -b master --recursive --depth 10
      )
    else()
      if(devMode)
        message("downloading latest Caliper updates")
        execute_process(COMMAND ${Git} pull origin master WORKING_DIRECTORY ${CALIPER_SRCDIR})
      endif(devMode)
    endif()

    option(CALIPER_OPTION_PREFIX ON)
    set(CALIPER_OPTION_PREFIX ON)

    option(CALIPER_WITH_MPI ON)
    set(CALIPER_WITH_MPI ON)

    option(ENABLE_TESTS "Enable Caliper Test" OFF )

    add_subdirectory(${CALIPER_SRCDIR})
    include_directories(${CMAKE_BINARY_DIR}/subprojects/caliper/include) # this is needed to find build-dir/include/SAMRAI/SAMRAI_config.h
    include_directories(${CALIPER_SRCDIR}/include)

  endif()

  add_definitions(-DPHARE_WITH_CALIPER=1)
  set (PHARE_BASE_LIBS ${PHARE_BASE_LIBS} caliper)
endif()
