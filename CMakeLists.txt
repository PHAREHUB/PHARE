cmake_minimum_required (VERSION 3.2)

project(phare VERSION 0.1)

set(CMAKE_CXX_STANDARD 14)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "" FORCE)
endif()


option(test "Build test with google test" ON)
option(coverage "Generate coverage" ON)
option(documentation "Add doxygen target to generate documentation" ON)

find_program(Git git)


find_package (PythonInterp  3.0 REQUIRED)

set(SUBPROJECTS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/subprojects)

set(CMAKE_EXPORT_COMPILE_COMMANDS 1)


if (test)
  set(GOOGLE_TEST_DIR ${CMAKE_CURRENT_SOURCE_DIR}/subprojects/googletest)

  if (NOT EXISTS ${GOOGLE_TEST_DIR})
    execute_process(
      COMMAND ${Git} clone https://github.com/google/googletest ${GOOGLE_TEST_DIR}
      )
  endif()
  add_subdirectory(subprojects/googletest)

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
                    '/usr/include/*' 'googletest/*' 'tests/*' '*/SAMRAI/*'         # that contain the pattern on this line
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

  # Test subdirectory goes here
endif()



if (documentation)

  find_package(Doxygen REQUIRED)

  set(DOXYGEN_GENERATE_HTML NO)
  set(DOXYGEN_GENERATE_LATEX NO)
  set(DOXYGEN_GENERATE_XML YES)

  doxygen_add_docs(
    doxygen
    ${CMAKE_SOURCE_DIR}/src  ${CMAKE_SOURCE_DIR}/tests
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    )

  set(doxy_dir ${CMAKE_CURRENT_BINARY_DIR}/xml)

  set(DOC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/subprojects/pharead)

  if (NOT EXISTS ${DOC_DIR})
    execute_process(
      COMMAND ${Git} clone https://hephaistos.lpp.polytechnique.fr/rhodecode/GIT_REPOSITORIES/LPP/phare/pharead
                     ${DOC_DIR}
      )
  endif()

  add_subdirectory(subprojects/pharead)


endif(documentation)

# Here goes src subdirectory

