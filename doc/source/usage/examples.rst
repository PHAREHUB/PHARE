
========
Examples
========

PHARE's functional test suite (``tests/functional/``) doubles as a set of
complete, runnable example scripts, each driving a real
``Simulation``/model/``Simulator`` setup for a specific physics problem.

Hybrid
------

* ``alfven_wave`` - propagation of a circularly polarized Alfvén wave.
* ``conservation`` - checks energy and momentum conservation over a run.
* ``dispersion`` - checks the Hybrid dispersion relation against theory, from noise and from prescribed modes.
* ``harris`` - magnetic reconnection in a Harris current sheet, 2D and 3D.
* ``ionIonBeam`` - ion-ion beam instability growth rate.
* ``shock`` - shock formation from a density/velocity jump.
* ``td`` / ``tdtagged`` - 1D tangential discontinuity, without and with AMR tagging.
* ``translation`` - advection of a uniform plasma / tangential discontinuity across the domain.

MHD
---

* ``mhd_alfven2d`` - 2D MHD Alfvén wave propagation.
* ``mhd_convergence`` / ``mhd_multidimensional_convergence`` - spatial convergence order of the MHD reconstruction schemes.
* ``mhd_dispersion`` - whistler wave dispersion relation.
* ``mhd_harris`` / ``mhd_harris_3d`` - MHD Harris current sheet reconnection, 2D and 3D.
* ``mhd_orszagtang`` / ``mhd_orszagtang_3d`` - Orszag-Tang vortex.
* ``mhd_rotor`` - MHD rotor problem.
* ``mhd_shock`` - MHD Riemann shock tube.
