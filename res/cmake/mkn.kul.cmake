

if(DEFINED MKN_KUL_ROOT)
  include_directories(${MKN_KUL_ROOT}/inc)
else()
  if (NOT EXISTS ${PHARE_PROJECT_DIR}/subprojects/mkn.kul)
    execute_process(
      COMMAND ${Git} clone https://github.com/mkn/mkn.kul ${PHARE_PROJECT_DIR}/subprojects/mkn.kul
    )
  endif()
  include_directories(
      ${PHARE_PROJECT_DIR}/subprojects/mkn.kul/inc
  )
endif()
