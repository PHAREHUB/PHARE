========
Examples
========

PHARE's functional tests (``tests/functional/``) are complete, runnable input
scripts, one per physics problem. Each one builds a ``Simulation``, a model and
some diagnostics, runs it, then plots or checks the result. They make good
starting points for your own runs. Each page below explains one example: the
physics, the setup, the exact command, and what the script checks.

How to run one
--------------

Set two shell variables: ``PHARE`` is your clone of PHARE and ``BUILD`` your
build folder (see :doc:`/getting_started/install` and
:doc:`/getting_started/build`). Then, in an **empty folder** (the scripts write
their output in the folder where you start them):

.. code-block:: bash

   export PHARE=/path/to/PHARE
   export BUILD=/path/to/PHARE/build
   mkdir -p my_run && cd my_run
   PYTHONPATH=$BUILD:$PHARE/pyphare:$PHARE \
       mpirun -n 4 python3 $PHARE/tests/functional/<name>/<script>.py

``$PHARE`` is on the path because some scripts import helpers from the
repository's ``tests`` package. The number of :term:`MPI <MPI>` processes is
free; every page says how many the test suite uses. To check that a script
only parses, without running it, set ``PHARE_DRY_RUN=1``.

Runtimes below were measured on a 24-core workstation with a release build,
on 4 processes, while other runs were using the machine. "Long run" means the
script was not run to the end here.

The examples
------------

.. list-table::
   :header-rows: 1
   :widths: 22 12 8 38 20

   * - Example
     - Model
     - Dim
     - Physics
     - Rough runtime
   * - :doc:`alfven_wave`
     - hybrid
     - 1D
     - circularly polarized Alfvén wave, phase speed
     - about 30 min (estimated)
   * - :doc:`ion_ion_beam`
     - hybrid
     - 1D
     - ion-ion beam instability, growth rate
     - about 20 min (estimated)
   * - :doc:`shock`
     - hybrid
     - 1D
     - smoothed magnetic jump, three interpolation orders
     - about 45 min (estimated)
   * - :doc:`conservation`
     - hybrid
     - 1D
     - kinetic energy drift in a uniform plasma
     - about 2 h (estimated)
   * - :doc:`tdtagged`
     - hybrid
     - 1D
     - advected tangential discontinuities, with and without AMR
     - 77 s
   * - :doc:`dispersion`
     - hybrid
     - 1D
     - whistler and Alfvén/ion-cyclotron dispersion relation
     - 117 s (one interpolation order)
   * - :doc:`harris`
     - hybrid
     - 2D, 3D
     - reconnection in Harris current sheets
     - 2D about 1 h (estimated); 3D long
   * - :doc:`mhd_alfven2d`
     - MHD
     - 2D
     - oblique Alfvén wave with AMR
     - 36 s
   * - :doc:`mhd_dispersion`
     - Hall MHD
     - 1D
     - whistler waves, two resolutions
     - Long run: not timed here
   * - :doc:`mhd_convergence`
     - MHD
     - 1D, 3D
     - convergence order of the reconstructions
     - Long run: not timed here
   * - :doc:`mhd_harris`
     - Hall MHD
     - 2D, 3D
     - reconnection in Harris current sheets, with boundaries
     - Long run: not timed here
   * - :doc:`mhd_orszagtang`
     - MHD
     - 2D, 3D
     - Orszag-Tang vortex
     - 2D: 102 s; 3D: long run
   * - :doc:`mhd_rotor`
     - MHD
     - 2D
     - rotor problem with AMR
     - 70 s
   * - :doc:`mhd_shock`
     - MHD
     - 1D
     - shock tube, periodic and with inflow boundary
     - 7 s and 15 s

.. toctree::
   :hidden:
   :maxdepth: 1

   alfven_wave
   ion_ion_beam
   shock
   conservation
   tdtagged
   dispersion
   harris
   mhd_alfven2d
   mhd_dispersion
   mhd_convergence
   mhd_harris
   mhd_orszagtang
   mhd_rotor
   mhd_shock
