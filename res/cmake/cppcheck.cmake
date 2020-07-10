
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
