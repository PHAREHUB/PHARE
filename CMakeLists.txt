cmake_minimum_required (VERSION 3.3)

project(phare VERSION 0.1)

set(CMAKE_CXX_STANDARD 17)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "" FORCE)
endif()

set(PHARE_PROJECT_DIR ${CMAKE_CURRENT_SOURCE_DIR})

option(test "Build test with google test" ON)
option(coverage "Generate coverage" ON)
option(documentation "Add doxygen target to generate documentation" ON)
option(cppcheck "Enable cppcheck xml report" ON)

option(asan "build with asan support" OFF)
option(ubsan "build with ubsan support" OFF)
option(msan "build with msan support" OFF)

find_program(Git git)

include(CheckCXXCompilerFlag)

find_package (PythonInterp  3.0 REQUIRED)

set(SUBPROJECTS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/subprojects)

set(CMAKE_EXPORT_COMPILE_COMMANDS 1)


#*******************************************************************************
#* Cppcheck option
#*******************************************************************************
if (cppcheck)
  find_program(Cppcheck cppcheck)
  if (NOT Cppcheck-NOTFOUND)

    add_custom_command(
      OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/cppcheck.xml
      COMMAND ${Cppcheck} --enable=all --std=c++11 --language=c++ --xml
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

  set(SAMRAI_SRCDIR ${CMAKE_CURRENT_SOURCE_DIR}/subprojects/samrai)
  set(SAMRAI_BIN ${CMAKE_CURRENT_BINARY_DIR}/subprojects/samrai)

  if (NOT EXISTS ${SAMRAI_SRCDIR})
    execute_process(
      COMMAND ${Git} clone https://github.com/monwarez/SAMRAI ${SAMRAI_SRCDIR} --branch add-install --recursive
      )
  endif()

  add_subdirectory(${SAMRAI_SRCDIR})


endif()



#*******************************************************************************
#* Test and coverage option
#*******************************************************************************
if (test)

  if (SAMRAI_FOUND)
      set(GOOGLE_TEST_DIR ${CMAKE_CURRENT_SOURCE_DIR}/subprojects/googletest)

      if (NOT EXISTS ${GOOGLE_TEST_DIR})
        execute_process(
          COMMAND ${Git} clone https://github.com/google/googletest ${GOOGLE_TEST_DIR}
          )
      endif()

   add_subdirectory(subprojects/googletest)
  endif()

  enable_testing()


  if (coverage)
    # For coverage  data we need to have zero optimization  and
    # to set special flags
    set(CMAKE_CXX_FLAGS "-g -O0 -Wall -fprofile-arcs -ftest-coverage")
    list(APPEND CMAKE_EXE_LINKER_FLAGS "-fprofile-arcs -ftest-coverage")

    # The followings is just to remove gcov extension, so that files will appear with their
    # actual extension
    set(CMAKE_CXX_OUTPUT_EXTENSION_REPLACE 1)

    # Here we made lcov aware of the entire source code
    # It will be invoked before test
    add_custom_command(
      OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/coverage-base.info                     # The data collected will be used as a base.
                                                                                # So that after test we can compute the
                                                                                # final coverage

      COMMAND lcov --capture --initial                                           # We capture initial data
                   --base-directory ${CMAKE_CURRENT_SOURCE_DIR}                  # The base source directory
                   --directory ${CMAKE_CURRENT_BINARY_DIR}                       # The binary directory
                   --output-file ${CMAKE_CURRENT_BINARY_DIR}/coverage-base.info  # Where we output the coverage
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      COMMENT "Pre-coverage report"
     )

     add_custom_command(
       OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/out/index.html

       COMMAND lcov --capture
                    --base-directory ${CMAKE_CURRENT_SOURCE_DIR}
                    --directory ${CMAKE_CURRENT_BINARY_DIR}
                    --output-file ${CMAKE_CURRENT_BINARY_DIR}/coverage-test.info

       COMMAND lcov
                    --add-tracefile ${CMAKE_CURRENT_BINARY_DIR}/coverage-base.info # Here we take coverage-base.info
                    --add-tracefile ${CMAKE_CURRENT_BINARY_DIR}/coverage-test.info # and coverage-test.info data
                    --output-file ${CMAKE_CURRENT_BINARY_DIR}/coverage.info        # and we merge them in coverage.info

       COMMAND lcov --remove ${CMAKE_CURRENT_BINARY_DIR}/coverage.info             # From coverage.info, we remove data
                    '/usr/include/*' '*/googletest/*' '*/tests/*' '*/SAMRAI/*'     # that contain the pattern on this line
                    --output-file ${CMAKE_CURRENT_BINARY_DIR}/filtered.info        # and put it in filtered.info

       COMMAND genhtml --ignore-errors source                                      # Since we remove data from coverage
                                                                                   # we need to tell genhtml to ignore
                                                                                   # error related of missing info from
                                                                                   # source files.
                       --prefix ${CMAKE_CURRENT_SOURCE_DIR}                        # Here we set the prefix of the repository
                       --keep-descriptions --frames --show-details --legend        # Some options for esthetic
                       --output-directory ${CMAKE_CURRENT_BINARY_DIR}/out          # Html output will be here
                       ${CMAKE_CURRENT_BINARY_DIR}/filtered.info                   # Tracefiles used (the filtered one)

        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating coverage repport"
       )
    add_custom_target(pre_coverage DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/coverage-base.info)
    add_custom_target(coverage_report DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/out/index.html)
  endif()


  add_subdirectory(tests/amr/data/particles)
  add_subdirectory(tests/amr/field_data)
  add_subdirectory(tests/amr/field_geometry)
  add_subdirectory(tests/amr/field_overlap)
  add_subdirectory(tests/amr/field_variable)
  add_subdirectory(tests/amr/tools/resources_manager)
  add_subdirectory(tests/core/data/ndarray)
  add_subdirectory(tests/core/data/field)
  add_subdirectory(tests/core/data/gridlayout)
  add_subdirectory(tests/core/data/vecfield)
  add_subdirectory(tests/core/data/particles)
  add_subdirectory(tests/core/data/ion_population)
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

endif()



#*******************************************************************************
#* Here the Core and AMR lib
#*******************************************************************************

add_subdirectory(src/core)
add_subdirectory(src/amr)

#*******************************************************************************
#* Documentation option
#*******************************************************************************
if (documentation)
  find_package(Doxygen REQUIRED)

  set(DOXYGEN_GENERATE_LATEX NO)
  set(DOXYGEN_GENERATE_XML YES)

  doxygen_add_docs(
    doxygen
    ${CMAKE_SOURCE_DIR}/src  ${CMAKE_SOURCE_DIR}/tests
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    )

  #set(doxy_dir ${CMAKE_CURRENT_BINARY_DIR}/xml)

  set(DOC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/subprojects/pharead)

  if (NOT EXISTS ${DOC_DIR})
    execute_process(
      COMMAND ${Git} clone https://hephaistos.lpp.polytechnique.fr/rhodecode/GIT_REPOSITORIES/LPP/phare/pharead
                     ${DOC_DIR}
      )
  endif()

  add_subdirectory(subprojects/pharead)


endif(documentation)


