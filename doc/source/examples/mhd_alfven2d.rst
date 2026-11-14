=======================
Alfvén wave in 2D (MHD)
=======================

**Goal:** propagate an Alfvén wave at 30 degrees to the grid axes with the
:term:`MHD` solver, on a grid with :term:`AMR`, and look at the result.

A uniform plasma (density 1, pressure 0.1) carries a mean field
:math:`(\cos\alpha, \sin\alpha)` with :math:`\alpha = 30^\circ`, and a circularly
polarized perturbation of amplitude 0.1 depending on
:math:`x\cos\alpha + y\sin\alpha`: the wave travels along the mean field. The
box size is :math:`1/\cos\alpha \times 1/\sin\alpha`, so that the wave is
periodic in both directions. The final time is 1, one wave period.

Setup
-----

.. list-table::
   :header-rows: 1

   * - Item
     - Value
   * - Dimension / domain
     - 2D, 100 x 100 cells, box 1.155 x 2, periodic
   * - Model
     - ideal MHD, :math:`\gamma = 5/3`, no Hall term, no resistivity or viscosity
   * - Scheme
     - linear reconstruction with the VanLeer limiter, Rusanov :term:`Riemann solver`, ``TVDRK2`` time stepping
   * - Time step
     - 0.002 (final time 1, 500 steps)
   * - :term:`AMR`
     - :term:`tagging`, up to 2 :term:`levels <level>`
   * - Output
     - B, density, velocity and pressure, at t = 0 and 1

.. literalinclude:: ../../../tests/functional/mhd_alfven2d/alfven2d.py
   :language: python
   :pyobject: config

How to run
----------

Work in an empty folder (see :doc:`/examples/index`). The test suite runs this
script on 4 :term:`MPI <MPI>` processes.

.. code-block:: bash

   mkdir -p mhd_alfven2d && cd mhd_alfven2d
   PYTHONPATH=$BUILD:$PHARE/pyphare:$PHARE \
       mpirun -n 4 python3 $PHARE/tests/functional/mhd_alfven2d/alfven2d.py

What it plots
-------------

There is no numerical assertion. At t = 0 and t = 1 the script saves PNG files
of :math:`\nabla\cdot B` (colour scale :math:`\pm 10^{-11}`), the MPI ranks,
:math:`v_x, v_y` and :math:`B_x, B_y`, with the patch outline, to
``phare_outputs/alfven2d_plots/<number of ranks>/``. The HDF5 files are
deleted at the end; only the plots stay.

Runtime
-------

36 s on 4 :term:`MPI <MPI>` processes (other runs were using the machine at
the same time).
