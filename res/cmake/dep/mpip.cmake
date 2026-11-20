
#  LLNL mpiP, lightweight MPI profiling via PMPI: https://github.com/LLNL/mpiP
#   1. configure cmake -DwithLLNLMpiP=ON
#      or -DMPIP_ROOT=/path/containing/libmpiP.so to use an existing build
#   2. needs libunwind and libbfd/libiberty headers, eg debian: libunwind-dev binutils-dev
#   3. python runs: mpirun -x LD_PRELOAD=${PHARE_MPIP_LIB} python3 sim.py
#      C++ executables/tests are linked with it directly
#   4. options via env, eg export MPIP="-k 4 -t 1.0", report written at MPI_Finalize
##

if(DEFINED MPIP_ROOT)
  set(withLLNLMpiP ON)
endif()

if (withLLNLMpiP)

  if(DEFINED MPIP_ROOT)
    find_library(PHARE_MPIP_LIB NAMES mpiP PATHS ${MPIP_ROOT} ${MPIP_ROOT}/lib NO_DEFAULT_PATH REQUIRED)
  else()

    if(NOT DEFINED PHARE_MPIP_VERSION)
      set(PHARE_MPIP_VERSION "v3.6.0")
    endif()

    set(MPIP_SRCDIR ${CMAKE_CURRENT_SOURCE_DIR}/subprojects/mpiP)
    set(PHARE_MPIP_LIB ${MPIP_SRCDIR}/libmpiP.so)

    if (NOT EXISTS ${MPIP_SRCDIR})
      message("cloning mpiP ${PHARE_MPIP_VERSION}")
      execute_process(COMMAND ${Git} clone -q --depth 1 -b ${PHARE_MPIP_VERSION}
                              https://github.com/LLNL/mpiP ${MPIP_SRCDIR}
                      COMMAND_ERROR_IS_FATAL ANY)
    endif()

    if (NOT EXISTS ${PHARE_MPIP_LIB})
      # partitioned: mpiP v3.6.0 MPI_Pready_list prototype conflicts with Open MPI 5.0 mpi.h
      message("building mpiP")
      execute_process(COMMAND ./configure --with-cc=${MPI_C_COMPILER} --disable-mpi-partitioned
                      WORKING_DIRECTORY ${MPIP_SRCDIR} OUTPUT_FILE ${MPIP_SRCDIR}/phare_configure.log
                      ERROR_FILE ${MPIP_SRCDIR}/phare_configure.log RESULT_VARIABLE MPIP_RESULT)
      if (NOT MPIP_RESULT EQUAL 0)
        message(FATAL_ERROR "mpiP configure failed, see ${MPIP_SRCDIR}/phare_configure.log"
                            " (libunwind/binutils dev packages installed?)")
      endif()
      execute_process(COMMAND make shared
                      WORKING_DIRECTORY ${MPIP_SRCDIR} OUTPUT_FILE ${MPIP_SRCDIR}/phare_make.log
                      ERROR_FILE ${MPIP_SRCDIR}/phare_make.log RESULT_VARIABLE MPIP_RESULT)
      if (NOT MPIP_RESULT EQUAL 0)
        message(FATAL_ERROR "mpiP build failed, see ${MPIP_SRCDIR}/phare_make.log")
      endif()
    endif()

  endif()

  message("mpiP for LD_PRELOAD: ${PHARE_MPIP_LIB}")
  set (PHARE_BASE_LIBS ${PHARE_BASE_LIBS} ${PHARE_MPIP_LIB})

endif(withLLNLMpiP)
