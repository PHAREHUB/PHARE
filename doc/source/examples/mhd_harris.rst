===========================================
Harris current sheets (Hall MHD, 2D and 3D)
===========================================

**Goal:** simulate magnetic reconnection in a pair of Harris current sheets
with the Hall :term:`MHD` solver, with periodic and with physical boundaries,
and in 3D.

Two current sheets sit at :math:`y = 0.25 L_y` and :math:`0.75 L_y`: :math:`B_x`
goes from -1 to 1 across the first and back across the second, over a width
0.5. The density is :math:`0.4 + \mathrm{sech}^2` of the distance to each sheet,
the plasma is at rest, and the pressure is :math:`p = 1 - (B_x^2 + B_y^2)/2`,
so that the total pressure is uniform. With boundaries the box is half as
tall, so it holds one sheet, at mid-height, between the two reflective walls.
In 3D the density is :math:`p/2` (uniform temperature) instead. A localized perturbation of
:math:`B_x, B_y` (amplitude 0.1) in the middle of the box triggers
reconnection.

Setup
-----

.. list-table::
   :header-rows: 1

   * - Item
     - periodic (``harris_periodic.py``)
     - with boundaries (``harris_with_boundaries.py``)
     - 3D (``harris_mhd_3d.py``)
   * - Cells / cell size
     - 160 x 80 / 0.375
     - 160 x 40 / 0.375
     - 150 x 240 x 50 / 0.4
   * - Boundaries
     - periodic
     - x: open; y: reflective
     - periodic
   * - Time step / final time
     - 0.005 / 50
     - 0.005 / 5
     - 0.002 / 0.1
   * - :term:`AMR`
     - :term:`tagging`, up to 2 :term:`levels <level>`
     - same
     - single level
   * - Scheme
     - ``WENOZ`` reconstruction, Rusanov :term:`Riemann solver`, ``TVDRK3``, Hall term on, :math:`\nu = 0.02`
     - same
     - same (the script also sets ``resistivity=0.001``, which has no effect in an MHD run; ``eta`` is 0)
   * - Output
     - B, density, velocity, pressure (6 times)
     - same
     - B and density at 0, 0.05, 0.1; VTK-HDF

.. literalinclude:: ../../../tests/functional/mhd_harris/common.py
   :language: python
   :pyobject: config

How to run
----------

Work in an empty folder (see :doc:`/examples/index`). The test suite runs all
three scripts on 10 :term:`MPI <MPI>` processes.

.. code-block:: bash

   mkdir -p mhd_harris && cd mhd_harris
   PYTHONPATH=$BUILD:$PHARE/pyphare:$PHARE \
       mpirun -n 10 python3 $PHARE/tests/functional/mhd_harris/harris_periodic.py
   PYTHONPATH=$BUILD:$PHARE/pyphare:$PHARE \
       mpirun -n 10 python3 $PHARE/tests/functional/mhd_harris/harris_with_boundaries.py
   PYTHONPATH=$BUILD:$PHARE/pyphare:$PHARE \
       mpirun -n 10 python3 $PHARE/tests/functional/mhd_harris_3d/harris_mhd_3d.py

What it checks and plots
------------------------

The 2D scripts save PNG files (:math:`\nabla\cdot B`, ranks, density, velocity,
:math:`B`, pressure, :math:`J_z`) for each output time to
``phare_outputs/<run>_plots/<number of ranks>/``. The periodic case has no
assertion. The case with boundaries raises an error if the output :math:`B`
contains any NaN, or if the wall-normal :math:`B_y` next to the two reflective
walls is larger than :math:`10^{-6}` at the final time. The 3D script only
writes its diagnostics.

Runtime
-------

Long run: not timed here (cluster recommended). The case with boundaries had
not finished after 180 s on 4 processes. All three scripts were checked in
dry-run mode.
