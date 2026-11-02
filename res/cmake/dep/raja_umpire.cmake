

if (withRaja)
  # to build see https://github.com/mkn-nix/llnl.raja mkn.sh
  find_package(raja REQUIRED
             NO_DEFAULT_PATH
             PATHS ${RAJA_DIR})

  set (PHARE_BASE_LIBS ${PHARE_BASE_LIBS} RAJA SAMRAI_tbox)
  set (PHARE_FLAGS ${PHARE_FLAGS} -DHAVE_RAJA=1)


endif(withRaja)


if (withUmpire)
  # to build see https://github.com/mkn-nix/llnl.umpire mkn.sh
  find_package(umpire REQUIRED
             NO_DEFAULT_PATH
             PATHS ${umpire_DIR})

  set (PHARE_BASE_LIBS ${PHARE_BASE_LIBS} umpire SAMRAI_tbox) # cudart
  set (PHARE_FLAGS ${PHARE_FLAGS} -DHAVE_UMPIRE=1)

endif(withUmpire)
