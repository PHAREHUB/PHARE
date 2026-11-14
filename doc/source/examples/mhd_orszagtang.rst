===================================
Orszag-Tang vortex (MHD, 2D and 3D)
===================================

**Goal:** run the Orszag-Tang vortex, a standard :term:`MHD` test of
interacting shocks and current sheets, in 2D and in a 3D version.

The initial state is smooth: uniform density :math:`25/(36\pi)` and pressure
:math:`5/(12\pi)`, a velocity :math:`v = (-\sin 2\pi y,\ \sin 2\pi x,\ 0)` and a
field :math:`B = B_0(-\sin 2\pi y,\ \sin 4\pi x,\ 0)` with
:math:`B_0 = 1/\sqrt{4\pi}`, in a periodic box of size 1. The 3D version adds a modulation :math:`1 + 0.2\sin 2\pi z` to
:math:`v_x, v_y` and a velocity :math:`v_z = 0.2\sin 2\pi z`.

Setup
-----

.. list-table::
   :header-rows: 1

   * - Item
     - 2D (``orszag_tang.py``)
     - 3D (``orszag_tang_3d.py``)
   * - Cells
     - 256 x 256
     - 64 x 64 x 64
   * - Model
     - ideal MHD, :math:`\gamma = 5/3`, no Hall term
     - same
   * - Scheme
     - linear reconstruction, VanLeer limiter, Rusanov :term:`Riemann solver`, ``TVDRK2``
     - same
   * - Time step
     - adaptive (:term:`CFL condition` number 0.8), final time 1
     - 0.00025, 100 steps (final time 0.025)
   * - :term:`AMR`
     - single :term:`level`
     - :term:`tagging`, up to 2 levels
   * - Output
     - B, density, velocity, pressure, 6 times from 0 to 1
     - same, 6 times

.. literalinclude:: ../../../tests/functional/mhd_orszagtang/orszag_tang.py
   :language: python
   :pyobject: config

How to run
----------

Work in an empty folder (see :doc:`/examples/index`). The test suite runs the 2D
script on 4 and the 3D script on 10 :term:`MPI <MPI>` processes.

.. code-block:: bash

   mkdir -p orszag_tang && cd orszag_tang
   PYTHONPATH=$BUILD:$PHARE/pyphare:$PHARE \
       mpirun -n 4 python3 $PHARE/tests/functional/mhd_orszagtang/orszag_tang.py
   PYTHONPATH=$BUILD:$PHARE/pyphare:$PHARE \
       mpirun -n 10 python3 $PHARE/tests/functional/mhd_orszagtang_3d/orszag_tang_3d.py

What it plots
-------------

There is no numerical assertion. The 2D script saves PNG files of
:math:`\nabla\cdot B`, the MPI ranks, density, :math:`v_x, v_y`,
:math:`B_x, B_y` and pressure for each time, to
``phare_outputs/orszag_tang_plots/<number of ranks>/``; the HDF5 files are
deleted at the end. The 3D script only runs. A third script in the 2D folder,
``orszag_tang_benchmark.py``, times the 2D case with a constant and with an
adaptive time step and prints the speed-up (no output files; checked in
dry-run mode only).

Runtime
-------

2D: 102 s on 4 :term:`MPI <MPI>` processes (other runs were using the machine
at the same time). 3D: long run, not timed here (cluster recommended): it was
at about 7 s per step on 4 processes.
