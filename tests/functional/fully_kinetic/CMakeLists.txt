cmake_minimum_required (VERSION 3.20.1)

project(test-fully-kinetic-waves)

if(NOT ${PHARE_PROJECT_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  file(GLOB PYFILES "*.py")
  file(COPY ${PYFILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif()

if(HighFive)
  ## These test use dump diagnostics so require HighFive!
  # Langmuir, light wave and two-stream, one test per interpolation order (~2.5 min each, Debug)
  foreach(order 1 2 3)
    phare_python3_exec(10 test-fully-kinetic-waves-o${order} waves.py ${CMAKE_CURRENT_BINARY_DIR} ${order})
    phare_test_labels(test-fully-kinetic-waves-o${order} physics)
  endforeach()
  # 2D: oblique Langmuir and light waves (both polarizations), order 1 only (cost);
  # orders 2 and 3: waves_2d.py 2 3, by hand
  phare_python3_exec(10 test-fully-kinetic-waves-2d-o1 waves_2d.py ${CMAKE_CURRENT_BINARY_DIR} 1)
  phare_test_labels(test-fully-kinetic-waves-2d-o1 physics)
  # 2D thermal plasma: the grid scale B_z growth is significant unfiltered and one current
  # filter pass cuts it (~46 s, Debug: two 16^2 runs); the longer runs of energy.py are by hand
  phare_python3_exec(10 test-fully-kinetic-energy-filter energy.py ${CMAKE_CURRENT_BINARY_DIR} --check)
  phare_test_labels(test-fully-kinetic-energy-filter physics)
  # self-force of one particle on one level against the Yee model, orders 1 to 3 (~30 s, Debug)
  phare_python3_exec(9 test-fully-kinetic-self-force self_force.py ${CMAKE_CURRENT_BINARY_DIR} --check)
  phare_test_labels(test-fully-kinetic-self-force physics)
  # the same at the refined runs' options (one filter pass, Marder s = 0.5), the reference the
  # refinement boundary gate subtracts (~12 s, Debug)
  phare_python3_exec(9 test-fully-kinetic-self-force-amr self_force.py ${CMAKE_CURRENT_BINARY_DIR} --check --filter-passes 1 --cleaning 0.5)
  phare_test_labels(test-fully-kinetic-self-force-amr physics)
  # T8: light pulse through a static 2:1 box, vacuum, reflection against the 1D Yee model,
  # orders 1 and 3, lambda = 32 and 16 fine cells (~40 s, Debug)
  phare_python3_exec(10 test-fully-kinetic-pulse-amr pulse.py ${CMAKE_CURRENT_BINARY_DIR} --check)
  phare_test_labels(test-fully-kinetic-pulse-amr physics)
  # the boundary node of a 1D fine level: no sawtooth between the band renewals (sawtooth fix),
  # order 1, 1600 particles per cell, 100 steps (165 s Debug, 1.9 GB of dumps)
  phare_python3_exec(10 test-fully-kinetic-sawtooth-gate sawtooth_gate.py ${CMAKE_CURRENT_BINARY_DIR})
  phare_test_labels(test-fully-kinetic-sawtooth-gate physics)
  # FK_rho dumped at fine time t(m+1) = the rho of the particles dumped at t(m) on every node of
  # the fine level of T2, boundary nodes included (1D, N = 0, s = 0.5, 6 steps)
  phare_python3_exec(9 test-fully-kinetic-fk-rho-pairing fk_rho_pairing.py ${CMAKE_CURRENT_BINARY_DIR})
  phare_test_labels(test-fully-kinetic-fk-rho-pairing physics)
  # BORDER-GAUSS: the Gauss residual of the edge nodes of a 2D fine level changes over a fine
  # substep as the deposit does, round-off, at N = 1 and 0 (no cleaning, 4 coarse steps, ~4 s Release)
  foreach(passes 1 0)
    phare_python3_exec(9 test-fully-kinetic-border-gauss-n${passes} border_gauss.py ${CMAKE_CURRENT_BINARY_DIR} --check --filter-passes ${passes})
    phare_test_labels(test-fully-kinetic-border-gauss-n${passes} physics)
  endforeach()
endif()
