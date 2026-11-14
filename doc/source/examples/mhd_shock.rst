==============
MHD shock tube
==============

**Goal:** run a 1D MHD Riemann problem, once in a periodic box and once with
an inflow boundary condition on the left.

The initial state has a jump in density, pressure and :math:`B_y`:

.. list-table::
   :header-rows: 1

   * - 
     - :math:`\rho`
     - :math:`p`
     - :math:`B_x`
     - :math:`B_y`
   * - left
     - 1
     - 1
     - 0.75
     - 1
   * - right
     - 0.125
     - 0.1
     - 0.75
     - -1

These are the states of the Brio-Wu shock tube (Brio and Wu 1988, which used
:math:`\gamma = 2`), here with :math:`\gamma = 5/3`.
The fluid is initially at rest in the periodic case, and moves at
:math:`1.05` times the largest fast magnetosonic speed of the two states in the other.

Setup
-----

.. list-table::
   :header-rows: 1

   * - Item
     - periodic (``shock_periodic.py``)
     - with boundaries (``shock_with_boundaries.py``)
   * - Cells
     - 800 of size 1 (box 800)
     - same
   * - Discontinuity at
     - :math:`x = 400`
     - :math:`x = 10`
   * - Boundaries
     - periodic
     - x lower: super-magnetofast inflow (left state imposed); x upper: open
   * - Time step / final time
     - 0.2 / 80 (400 steps)
     - 0.0667 / 66.7 (1000 steps)
   * - Scheme
     - ideal MHD, linear reconstruction, VanLeer limiter, Rusanov :term:`Riemann solver`, ``TVDRK2``
     - same
   * - :term:`AMR`
     - no, a single :term:`level`
     - no
   * - Output
     - B, density, velocity, pressure at the end
     - the same, every 100 steps

.. literalinclude:: ../../../tests/functional/mhd_shock/common.py
   :language: python
   :pyobject: config

How to run
----------

Work in an empty folder (see :doc:`/examples/index`). The test suite runs both
scripts on 10 :term:`MPI <MPI>` processes.

.. code-block:: bash

   mkdir -p mhd_shock && cd mhd_shock
   PYTHONPATH=$BUILD:$PHARE/pyphare:$PHARE \
       mpirun -n 4 python3 $PHARE/tests/functional/mhd_shock/shock_periodic.py
   PYTHONPATH=$BUILD:$PHARE/pyphare:$PHARE \
       mpirun -n 4 python3 $PHARE/tests/functional/mhd_shock/shock_with_boundaries.py

What it checks and plots
------------------------

Both scripts save PNG files of density, velocity, :math:`B_y` and pressure to
``phare_outputs/<run>_plots/<number of ranks>/``. The periodic one has no
assertion. The inflow one asserts that the density is finite everywhere on
level 0 and that the first interior cell next to the inflow boundary keeps
the left density within 5%.

Runtime
-------

Periodic: 7 s; with boundaries: 15 s, both on 4 :term:`MPI <MPI>` processes.
