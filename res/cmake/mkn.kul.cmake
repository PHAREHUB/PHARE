
if (NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/subprojects/mkn.kul)

  execute_process(
    COMMAND ${Git} clone https://github.com/mkn/mkn.kul ${CMAKE_CURRENT_SOURCE_DIR}/subprojects/mkn.kul}
  )

endif()

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/subprojects/mkn.kul/inc
)
