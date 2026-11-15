==============================
MHD validation: what is tested
==============================

**Goal:** know, for each MHD test of the repository, what it really checks, whether it runs
by default, what was measured at the commit documented here, and what is not tested.

The set-ups themselves are on the :doc:`examples </examples/index>` pages
(:doc:`shock </examples/mhd_shock>`, :doc:`Orszag–Tang </examples/mhd_orszagtang>`,
:doc:`rotor </examples/mhd_rotor>`, :doc:`Harris </examples/mhd_harris>`,
:doc:`Alfvén wave </examples/mhd_alfven2d>`, :doc:`convergence </examples/mhd_convergence>`,
:doc:`dispersion </examples/mhd_dispersion>`). The numerics are on
:doc:`/numerics/mhd_finite_volume`.

.. contents:: On this page
   :local:
   :depth: 1

What runs by default
--------------------

Each test is registered at an *execution level*. A default build runs the levels 1 to 10
only (:src:`res/cmake/options.cmake#L84-L90`); the functional MHD tests are at levels 11,
21, 101 and 121, so none of them is in a default ``ctest``. To run them, configure with a
higher ``PHARE_EXEC_LEVEL_MAX`` and with ``-DtestMPI=ON`` (off by default,
:src:`res/cmake/options.cmake#L15`); they also need HighFive for the diagnostics. The test
names below are those of a ``testMPI=ON`` build; without it the MPI suffix is absent.

In a default ``ctest``, the MHD checks are: the initial magnetic field in 1D and the
agreement of shared faces between patches, and, in the default serial configuration, the
Python interface without hybrid keys and the time interpolation of regridded ghosts. None of them advances an MHD simulation and checks
its accuracy.

The repository has a CI script that configures levels 11 to 100 with ``PHARE_DRY_RUN=1``
(:src:`res/ci/jobs/pull_request/dry_run_functional_tests.sh#L11-L12`). In a dry run the
simulation is built and validated but neither initialised nor advanced
(:src:`pyphare/pyphare/pharein/simulation.py#L1145`).

Tests that check something
--------------------------

* **1D convergence** (``py3_test-mhd_convergence_mpi_n_4``, level 21, 4 ranks). Runs each
  of the five reconstructions on a small Alfvén wave for several grids and asserts that the
  fitted slope of the error is within 15 % of 1, 2, 3, 5, 5
  (:src:`tests/functional/mhd_convergence/convergence.py#L170-L174`).
* **Alfvén and whistler waves with AMR** (``py3_amr-convergence-alfven2d_mpi_n_2`` and
  ``py3_amr-convergence-whistler2d_mpi_n_2``, level 21, 2D, ``WENOZ``, ``SSPRK4_5``). The
  measured order of the error must lie in the band 1.70 to 2.25, and the ratio of the AMR
  error to the uniform-grid error must not drift by more than 0.05 over a sweep
  (:src:`tests/simulator/amr_convergence/test_alfven2d_amr.py#L49`,
  :src:`tests/simulator/amr_convergence/test_alfven2d_amr.py#L58`). This is the only MHD
  test that measures the order of the scheme in more than one dimension.
* **Shock tube with boundaries** (``py3_mhd_shock_with_boundaries_mpi_n_10``, level 11,
  10 ranks). The density is finite and the first interior cell stays within 5 % of the
  left state (:src:`tests/functional/mhd_shock/shock_with_boundaries.py#L66`,
  :src:`tests/functional/mhd_shock/shock_with_boundaries.py#L76`).
* **Harris sheet with boundaries** (``py3_mhd_harris_with_boundaries_mpi_n_10``, level 11).
  No NaN in B, and :math:`|B_y|` at the walls below :math:`10^{-6}`
  (:src:`tests/functional/mhd_harris/harris_with_boundaries.py#L74-L87`).
* **Initial B in 1D** (level 9): B equals the user function to :math:`10^{-15}`; this is the
  initialisation only.
* **Shared faces** (level 9, 2D, 3 levels, 2 steps): B is identical, bit for bit, on faces
  shared by two patches.
* **Regridded ghosts**: ghost values at regrid time are linear in time to :math:`10^{-12}`.

Tests that only run or plot
---------------------------

These tests have no assertion that can fail on the numerics; they can fail only by
crashing or hanging. They are useful as demonstrations and to catch crashes.

* the Alfvén wave in 2D, the Orszag–Tang vortex in 2D (levels 21 and 101 for the
  benchmark) and 3D, the rotor, the whistler wave in 1D, the periodic Harris sheet and the
  periodic shock tube (plots or diagnostics only);
* the 3D Harris sheet: see the note below.

Measured at this commit
-----------------------

These are *documentation runs*, made with the scripts below on the code commit named in the header of each data file. They
are not part of the test suite.

**Convergence of each reconstruction** (:doc:`/numerics/mhd_reconstruction`, with the
figure): the 1D wave of the test, amplitude :math:`10^{-6}`, ``SSPRK4_5``, time step
:math:`5\times10^{-4}`, ``Rusanov``. The local orders are 0.80 to 0.97 (``Constant``),
1.84 to 2.03 (``Linear`` + ``VanLeer``), 1.70 to 1.93 (``Linear`` + ``MinMod``), 2.88 to
2.99 (``WENO3``), 4.82 to 4.86 (``WENOZ``) and 4.63 to 4.85 (``MP5``). They show the order
of each reconstruction on a small linear wave; they are not the order of the scheme, which
is second order at best.

**Order in time** (:doc:`/numerics/mhd_time_integration`): 1 for ``Euler``, 2 for
``TVDRK2``, 3 for ``TVDRK3``; the order 4 of ``SSPRK4_5`` is checked on its coefficients.

**Diffusion of the Riemann solvers** (:doc:`/numerics/mhd_riemann_solvers`): width of the
contact discontinuity of the shock tube, 45, 35 and 21 cells for ``Rusanov``, ``HLL`` and
``HLLD`` with ``Constant`` reconstruction, 7 cells with ``Linear`` + ``VanLeer`` + ``Rusanov``.

Not tested
----------

* The result of ``HLL``, ``HLLD`` and of the ``MinMod`` limiter. They are compiled
  (:src:`res/sim/all.txt#L63-L66`, "coverage-only") and no test checks their output; all
  the test decks use ``Rusanov``.
* There is no unit test of a reconstruction, a limiter, a Riemann solver or a time
  integrator.
* The order of the scheme in 3D.
* The MHD advance of the simulator-level field tests (the MHD variant of the
  advance tests is never selected).

Known issues
------------

.. note:: **Single-precision output limits the 1D convergence test.**
   The diagnostics are written in single precision, so errors under about
   :math:`10^{-13}` cannot be measured. For ``WENOZ`` and ``MP5`` the error reaches
   exactly zero on the finest grids of the test, the fit of the slope is then not a
   number, and the check for ``WENOZ`` fails ("got nan"). The documentation runs above
   use coarser grids and ten periods for these two reconstructions to stay above this
   floor.

.. note:: **The 3D convergence test does not run a simulation.** The script of
   ``mhd_multidimensional_convergence`` has its ``run()`` call commented out
   (:src:`tests/functional/mhd_multidimensional_convergence/multidimensional_convergence.py#L254`);
   it reads output that must already be there. It is registered at level 121.

.. note:: **The positivity check of the 3D Harris sheet is not called.** The assertion that
   the temperature is positive sits in a helper function that nothing calls
   (:src:`tests/functional/mhd_harris_3d/harris_mhd_3d.py#L111`).

Scripts and data
----------------

The data files are in ``doc/source/validation/data/``; each starts with a header giving the
commit, build, command and date. One script runs everything and rewrites them:

.. code-block:: bash

   python3 doc/source/validation/scripts/run_doc_runs.py --build <PHARE build dir>

It takes about 40 seconds on a workstation (33 serial runs, four at a time by default,
``--jobs`` changes it; about 20 seconds with eight). With
``--check`` it reruns into a temporary folder and compares with the committed files, with a
relative tolerance for each quantity; errors at the precision floor of the single-precision
output are not compared, so a build that writes double precision, or another compiler, can
be checked. The exit status is non-zero on any drift. ``--list-modules`` prints the PHARE
modules the runs need: in 1D, with ``Rusanov`` unless stated, ``Constant``, ``Linear`` +
``VanLeer``, ``Linear`` + ``MinMod``, ``WENO3``, ``WENOZ``, ``MP5``, and ``Constant`` with
``HLL`` and ``HLLD``; in 2D, ``Linear`` + ``VanLeer`` + ``Rusanov`` (the divergence runs).

.. literalinclude:: scripts/run_doc_runs.py
   :language: python
   :lines: 1-20
   :caption: run_doc_runs.py: runs everything, rewrites or checks the data files (start of the file)

.. literalinclude:: scripts/mhd_convergence_scan.py
   :language: python
   :caption: mhd_convergence_scan.py: convergence of one reconstruction

Other scripts: :download:`run_doc_runs.py <scripts/run_doc_runs.py>`,
:download:`mhd_shock_solvers.py <scripts/mhd_shock_solvers.py>`,
:download:`mhd_shock_width.py <scripts/mhd_shock_width.py>`,
:download:`mhd_time_order_ratios.py <scripts/mhd_time_order_ratios.py>`,
:download:`ssprk_check.py <scripts/ssprk_check.py>`.
