function(phare_cgroup_limit_mb file out_var)
  set(${out_var} "" PARENT_SCOPE)
  if(NOT EXISTS "${file}")
    return()
  endif()
  file(READ "${file}" content)
  string(STRIP "${content}" content)
  if("${content}" MATCHES "^[0-9]+$")
    math(EXPR limit_mb "${content} / 1048576")
    set(${out_var} "${limit_mb}" PARENT_SCOPE)
  endif()
endfunction()

function(phare_module_compile_jobs total_mb cgroup_mb per_module_mb out_var)
  set(memory_mb ${total_mb})
  if(NOT "${cgroup_mb}" STREQUAL "" AND cgroup_mb LESS memory_mb)
    set(memory_mb ${cgroup_mb})
  endif()
  math(EXPR jobs "${memory_mb} / ${per_module_mb}")
  if(jobs LESS 1)
    set(jobs 1)
  endif()
  set(${out_var} ${jobs} PARENT_SCOPE)
endfunction()

if(CMAKE_SCRIPT_MODE_FILE)
  return()
endif()

cmake_host_system_information(RESULT PHARE_HOST_MEMORY_MB QUERY TOTAL_PHYSICAL_MEMORY)
phare_cgroup_limit_mb("/sys/fs/cgroup/memory.max" PHARE_CGROUP_MEMORY_MB)
if("${PHARE_CGROUP_MEMORY_MB}" STREQUAL "")
  phare_cgroup_limit_mb("/sys/fs/cgroup/memory/memory.limit_in_bytes" PHARE_CGROUP_MEMORY_MB)
endif()

set(PHARE_MODULE_COMPILE_MEMORY_MB "2048" CACHE STRING "Peak memory in MiB of one simulator module compile")
phare_module_compile_jobs(
  ${PHARE_HOST_MEMORY_MB} "${PHARE_CGROUP_MEMORY_MB}" ${PHARE_MODULE_COMPILE_MEMORY_MB} PHARE_DEFAULT_MODULE_COMPILE_JOBS)
set(PHARE_MODULE_COMPILE_JOBS "${PHARE_DEFAULT_MODULE_COMPILE_JOBS}" CACHE STRING "Maximum concurrent simulator module compiles, 0 disables the limit")

message(STATUS "PHARE module compiles: at most ${PHARE_MODULE_COMPILE_JOBS} at once (${PHARE_HOST_MEMORY_MB} MiB, ${PHARE_MODULE_COMPILE_MEMORY_MB} MiB each)")

if(CMAKE_GENERATOR MATCHES "Ninja" AND PHARE_MODULE_COMPILE_JOBS GREATER 0)
  set_property(GLOBAL APPEND PROPERTY JOB_POOLS phare_modules=${PHARE_MODULE_COMPILE_JOBS})
endif()

function(phare_limit_module_compile target)
  if(NOT PHARE_MODULE_COMPILE_JOBS GREATER 0)
    return()
  endif()
  if(CMAKE_GENERATOR MATCHES "Ninja")
    set_property(TARGET ${target} PROPERTY JOB_POOL_COMPILE phare_modules)
  elseif(CMAKE_GENERATOR MATCHES "Makefiles")
    set(launcher
      ${Python_EXECUTABLE}
      ${PHARE_PROJECT_DIR}/res/cmake/compile_slot.py
      ${PHARE_MODULE_COMPILE_JOBS}
      ${CMAKE_BINARY_DIR}/.phare_compile_slots)
    if(CMAKE_CXX_COMPILER_LAUNCHER)
      list(APPEND launcher ${CMAKE_CXX_COMPILER_LAUNCHER})
    endif()
    set_property(TARGET ${target} PROPERTY CXX_COMPILER_LAUNCHER "${launcher}")
  endif()
endfunction()
