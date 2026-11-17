=========
MHD rotor
=========

**Goal:** run the MHD rotor problem, a dense disc spinning inside a lighter
medium threaded by a uniform magnetic field, with :term:`MHD` and :term:`AMR`.

The box is :math:`1 \times 1`, the field :math:`B = (B_0, 0, 0)` with
:math:`B_0 = 5/\sqrt{4\pi}`, the pressure 1. Inside a radius 0.1 around the
centre the density is 10 and the fluid rotates rigidly (speed 2 at the edge);
between radius 0.1 and 0.115 the density and the speed decrease linearly to
the background (density 1, at rest).

Setup
-----

.. list-table::
   :header-rows: 1

   * - Item
     - Value
   * - Dimension / domain
     - 2D, 100 x 100 cells of size 0.01, periodic
   * - Model
     - ideal MHD, :math:`\gamma = 5/3`, no Hall term
   * - Scheme
     - linear reconstruction, VanLeer limiter, Rusanov :term:`Riemann solver`, ``TVDRK2``
   * - Time step
     - adaptive (:term:`CFL condition` number 0.8), final time 0.15
   * - :term:`AMR`
     - :term:`tagging`, up to 2 :term:`levels <level>`
   * - Output
     - B, density, velocity, pressure, 6 times from 0 to 0.15

.. literalinclude:: ../../../tests/functional/mhd_rotor/rotor.py
   :language: python
   :pyobject: config

How to run
----------

Work in an empty folder (see :doc:`/examples/index`). The test suite runs this
script on 2 :term:`MPI <MPI>` processes.

.. code-block:: bash

   mkdir -p mhd_rotor && cd mhd_rotor
   PYTHONPATH=$BUILD:$PHARE/pyphare:$PHARE \
       mpirun -n 4 python3 $PHARE/tests/functional/mhd_rotor/rotor.py

What it plots
-------------

There is no numerical assertion. For each time the script saves PNG files of
:math:`\nabla\cdot B`, the MPI ranks, density, :math:`v_x, v_y`,
:math:`B_x, B_y` and pressure, with patch outlines, to
``phare_outputs/rotor_plots/<number of ranks>/``. The HDF5 files are deleted
at the end.

Runtime
-------

70 s on 4 :term:`MPI <MPI>` processes (other runs were using the machine at the
same time).
