cmake_minimum_required (VERSION 3.3)

project(PHARE VERSION 0.1)

set(CMAKE_CXX_STANDARD 17)
if (POLICY CMP0074) # hides warning about ${PACKAGE}_ROOT variables
  cmake_policy(SET CMP0074 NEW)
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "" FORCE)
endif()

set(PHARE_PROJECT_DIR ${CMAKE_CURRENT_SOURCE_DIR})

option(HighFive "Build with highfive usage" ON)
option(test "Build test with google test" ON)
option(coverage "Generate coverage" ON)
option(documentation "Add doxygen target to generate documentation" ON)
option(cppcheck "Enable cppcheck xml report" ON)
option(asan "build with asan support" OFF)
option(ubsan "build with ubsan support" OFF)
option(msan "build with msan support" OFF)

find_program(Git git)

include(CheckCXXCompilerFlag)
include(ExternalProject)

find_package (PythonInterp  3.0 REQUIRED)

set(SUBPROJECTS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/subprojects)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

# Pybind errors with clang, it is default in GCC
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set (CMAKE_CXX_FLAGS -fsized-deallocation)
endif()
if(MSVC)
  set (PHARE_WERROR_FLAGS /W4 /WX)
else()
  set (PHARE_WERROR_FLAGS -Wall -Wextra -pedantic -Werror)
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set (PHARE_WERROR_FLAGS ${PHARE_WERROR_FLAGS}
                        -Wno-unused-variable -Wno-unused-parameter
                        -Wno-gnu-zero-variadic-macro-arguments)
else()
  set (PHARE_WERROR_FLAGS ${PHARE_WERROR_FLAGS}
                          -Wno-unused-variable -Wno-unused-parameter
                          -Wno-unused-but-set-variable -Wno-unused-but-set-parameter)
endif()
endif()

#*******************************************************************************
#* Cppcheck option
#*******************************************************************************
if (cppcheck)
  find_program(Cppcheck cppcheck)
  if (NOT Cppcheck-NOTFOUND)

    add_custom_command(
      OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/cppcheck.xml
      COMMAND ${Cppcheck} --enable=all --std=c++11 --language=c++ --xml -i${CMAKE_CURRENT_LIST_DIR}/subprojects
      --project=${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json 2> ${CMAKE_CURRENT_BINARY_DIR}/cppcheck.xml
      )

    add_custom_target(cppcheck-xml DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/cppcheck.xml)
    find_program(Cppcheck-html cppcheck-htmlreport)
    if (NOT Cppcheck-html-NOTFOUND)

      add_custom_command(
        OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/cppcheckHtml/index.html
        DEPENDS cppcheck-xml
        COMMAND ${Cppcheck-html} --file=${CMAKE_CURRENT_BINARY_DIR}/cppcheck.xml
                                 --report-dir=${CMAKE_CURRENT_BINARY_DIR}/cppcheckHtml
        )

      add_custom_target(cppcheck-html DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/cppcheckHtml/index.html)
    endif()

  endif()
endif()

#*******************************************************************************
#* SAMRAI build option
#*******************************************************************************
find_package(SAMRAI CONFIG QUIET)
if (NOT SAMRAI_FOUND)

  if(DEFINED SAMRAI_ROOT)
    find_package(SAMRAI PATHS ${SAMRAI_ROOT} REQUIRED)
  else()
    set(SAMRAI_SRCDIR ${CMAKE_CURRENT_SOURCE_DIR}/subprojects/samrai)
    set(SAMRAI_BIN ${CMAKE_CURRENT_BINARY_DIR}/subprojects/samrai)

    if (NOT EXISTS ${SAMRAI_SRCDIR})
      execute_process(
        COMMAND ${Git} clone https://github.com/LLNL/SAMRAI ${SAMRAI_SRCDIR} -b master --recursive --depth 10
        )
    endif()

    option(ENABLE_TESTS "Enable Samrai Test" OFF ) # disable SAMRAI Test so that we can use the googletest pulled after

    add_subdirectory(${SAMRAI_SRCDIR})
  endif()

endif()

if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/build/include/SAMRAI")
include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/build/include
)
endif()

if(HighFive)
  include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/subprojects/highfive/include
  )

  if(DEFINED HDF5_ROOT)
    find_package(HDF5 PATHS ${HDF5_ROOT} REQUIRED)
  else()
    find_package(HDF5 REQUIRED)
  endif()

  message("HighFive enabled - checking HDF5")
  if(${HDF5_IS_PARALLEL})
      message("HDF5 PARALLEL detected")
      add_definitions(-DH5_HAVE_PARALLEL)
  else()
      message(WARNING "HDF5 NOT PARALLEL")
  endif()

  message(STATUS "HDF5_LIBRARIES " ${HDF5_LIBRARIES})
  message(STATUS "HDF5_INCLUDE_DIRS " ${HDF5_INCLUDE_DIRS})
  message(STATUS "HDF5_LIBRARY_PATH " ${HDF5_LIBRARY_PATH})
endif()

#*******************************************************************************
#* Pybind11 build option
#*******************************************************************************
find_package(pybind11 CONFIG QUIET)

if (NOT pybind11_FOUND)

    message("downloading subproject pybind11")
    set(PYBIND11_SRCDIR ${CMAKE_CURRENT_SOURCE_DIR}/subprojects/pybind11)

    if (NOT EXISTS ${PYBIND11_SRCDIR})
      execute_process(
        COMMAND ${Git} clone https://github.com/pybind/pybind11 ${PYBIND11_SRCDIR}
        )
    endif()

    add_subdirectory(${PYBIND11_SRCDIR})

endif()

#*******************************************************************************
#* Test and coverage option
#*******************************************************************************
if (test)

  set(GOOGLE_TEST_DIR ${CMAKE_CURRENT_SOURCE_DIR}/subprojects/googletest)

  if (NOT EXISTS ${GOOGLE_TEST_DIR})
     execute_process(
     COMMAND ${Git} clone https://github.com/google/googletest ${GOOGLE_TEST_DIR}
     )
  endif()

  add_subdirectory(subprojects/googletest)

  enable_testing()

  IF(coverage)
      set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -g -O0 -Wall -W -Wshadow -Wunused-variable -Wunused-parameter -Wunused-function -Wunused -Wno-system-headers -Wno-deprecated -Woverloaded-virtual -Wwrite-strings -fprofile-arcs -ftest-coverage")
      set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -g -O0 -Wall -W -Wshadow -Wunused-variable \
          -Wunused-parameter -Wunused-function -Wunused -Wno-system-headers \
          -Wno-deprecated -Woverloaded-virtual -Wwrite-strings -fprofile-arcs -ftest-coverage")
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fprofile-arcs -ftest-coverage")

      add_custom_target(build-time-make-directory ALL
          COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/coverage)

      add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/coverage/index.html
          COMMAND gcovr --exclude='.*subprojects.*' --exclude='.*tests.*' --exclude='/usr/include/.*' --object-directory ${CMAKE_BINARY_DIR}  -r ${CMAKE_SOURCE_DIR} --html  --html-details -o ${CMAKE_CURRENT_BINARY_DIR}/coverage/index.html
          )

    add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/coverage/coverage.xml
          COMMAND gcovr --exclude='.*subprojects.*' --exclude='.*tests.*' --exclude='/usr/include/.*' --object-directory ${CMAKE_BINARY_DIR}  -r ${CMAKE_SOURCE_DIR} --xml -o ${CMAKE_CURRENT_BINARY_DIR}/coverage/coverage.xml
          )

      add_custom_target(gcovr
          DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/coverage/index.html gcovr ${CMAKE_CURRENT_BINARY_DIR}/coverage/coverage.xml
      )

      if(APPLE)
          set(OPPEN_CMD open)
      elseif(UNIX)
          set(OPPEN_CMD xdg-open)
      endif(APPLE)

      add_custom_target(show_coverage
          COMMAND ${OPPEN_CMD} ${CMAKE_CURRENT_BINARY_DIR}/coverage/index.html
          DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/coverage/index.html gcovr
          )
  ENDIF(coverage)

  include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/tests
  )

  add_subdirectory(tests/amr/data/particles)
  add_subdirectory(tests/amr/data/field/coarsening)
  add_subdirectory(tests/amr/data/field/copy_pack)
  add_subdirectory(tests/amr/data/field/geometry)
  add_subdirectory(tests/amr/data/field/overlap)
  add_subdirectory(tests/amr/data/field/refine)
  add_subdirectory(tests/amr/data/field/variable)
  add_subdirectory(tests/amr/data/field/time_interpolate)
  add_subdirectory(tests/amr/resources_manager/)
  add_subdirectory(tests/amr/messengers)
  add_subdirectory(tests/amr/models)
  add_subdirectory(tests/amr/multiphysics_integrator)
  add_subdirectory(tests/core/data/ndarray)
  add_subdirectory(tests/core/data/field)
  add_subdirectory(tests/initializer)
  add_subdirectory(tests/core/data/gridlayout)
  add_subdirectory(tests/core/data/vecfield)
  add_subdirectory(tests/core/data/particles)
  add_subdirectory(tests/core/data/ions)
  add_subdirectory(tests/core/data/ion_population)
  add_subdirectory(tests/core/data/maxwellian_particle_initializer)
  add_subdirectory(tests/core/data/particle_initializer)
  add_subdirectory(tests/core/utilities/box)
  add_subdirectory(tests/core/utilities/particle_selector)
  add_subdirectory(tests/core/utilities/partitionner)
  add_subdirectory(tests/core/utilities/range)
  add_subdirectory(tests/core/utilities/index)
  add_subdirectory(tests/core/numerics/boundary_condition)
  add_subdirectory(tests/core/numerics/interpolator)
  add_subdirectory(tests/core/numerics/pusher)
  add_subdirectory(tests/core/numerics/ampere)
  add_subdirectory(tests/core/numerics/faraday)
  add_subdirectory(tests/core/numerics/ohm)
  add_subdirectory(tests/diagnostic)

endif()






#*******************************************************************************
#* Build the different libs and executables
#*******************************************************************************

add_subdirectory(src/core)
add_subdirectory(src/solver)
add_subdirectory(src/diagnostic)
add_subdirectory(src/initializer)
add_subdirectory(src/phare)
add_subdirectory(src/pyphare)
add_subdirectory(src/amr)





#*******************************************************************************
#* Documentation option
#*******************************************************************************
if (documentation)
    find_package(Doxygen)
    option(BUILD_DOCUMENTATION "Create and install the HTML based API documentation (requires Doxygen)" ${DOXYGEN_FOUND})

    if(BUILD_DOCUMENTATION)
        if(NOT DOXYGEN_FOUND)
            message(FATAL_ERROR "Doxygen is needed to build the documentation.")
        endif()

        set(doxyfile_in ${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile.in)
        set(doxyfile ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

        # The project version number.
        set(VERSION_MAJOR   1  )#CACHE STRING "Project major version number.")
        set(VERSION_MINOR   0  )#CACHE STRING "Project minor version number.")
        set(VERSION_PATCH   0  )#CACHE STRING "Project patch version number.")
        set(doxy_main_page "doc/phare.md")
        mark_as_advanced(VERSION_MAJOR VERSION_MINOR VERSION_PATCH)

        configure_file(${doxyfile_in} ${doxyfile} @ONLY)

        add_custom_target(doxygen
            COMMAND ${DOXYGEN_EXECUTABLE} ${doxyfile}  ${CMAKE_SOURCE_DIR}/doc ${CMAKE_SOURCE_DIR}/src  ${CMAKE_SOURCE_DIR}/tests
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Generating API documentation with Doxygen"
            VERBATIM)

        install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html DESTINATION share/doc)
    endif()



  set(DOC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/subprojects/pharead)

  if (NOT EXISTS ${DOC_DIR})
    execute_process(
      COMMAND ${Git} clone https://hephaistos.lpp.polytechnique.fr/rhodecode/GIT_REPOSITORIES/LPP/phare/pharead
                     ${DOC_DIR}
      )
  endif()

  add_subdirectory(subprojects/pharead)


endif(documentation)


#*******************************************************************************
#* Install PHAREIN
#*******************************************************************************
set(PHAREIN_DIR "${PROJECT_SOURCE_DIR}/external/pharein")


