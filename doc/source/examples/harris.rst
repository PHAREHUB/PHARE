=======================================================
Harris current sheets: reconnection (hybrid, 2D and 3D)
=======================================================

**Goal:** simulate magnetic reconnection in a pair of Harris current sheets
with the :term:`hybrid PIC` model, in 2D and in 3D.

Two current sheets sit at :math:`y = 0.3 L_y` and :math:`0.7 L_y`, where the
:math:`B_x` field changes sign (from -1 to 1 across the first, back across the
second). The density is :math:`0.4 + \mathrm{sech}^2` of the distance to each
sheet, with a width of 0.5. The temperature is set so that the total pressure
is uniform: :math:`T = (0.7 - B^2/2)/n`. Two small localized perturbations of
:math:`B_x, B_y` (amplitude 0.1, width 1), one at the middle of each sheet,
trigger reconnection. Electrons are isothermal with :math:`T_e = 0`. The script sets
the initial state and runs; it does not compare with a formula.

Setup
-----

.. list-table::
   :header-rows: 1

   * - Item
     - 2D (``harris_2d.py``)
     - 3D (``harris_3d.py``)
   * - Cells / cell size
     - 200 x 100 / 0.4 (box 80 x 40)
     - 40 x 80 x 20 / 0.4 (box 16 x 32 x 8)
   * - Model, particles
     - hybrid, protons, 100 per cell
     - hybrid, protons, 100 per cell
   * - Time step / final time
     - 0.005 / 50 (10 000 steps)
     - 0.001 / 1 (1000 steps)
   * - Resistivity
     - 0.001, hyper-resistivity 0.002
     - 0.001, hyper-resistivity 0.008 (spatial)
   * - :term:`AMR`
     - :term:`tagging`, up to 2 :term:`levels <level>`
     - :term:`tagging` is on but ``max_nbr_levels=1``: a single level
   * - Output
     - E, B, mass density, bulk velocity, proton density and pressure tensor, at 0, 25 and 50; HDF5 (``phareh5``)
     - B, E, bulk velocity, proton density, at the final time; VTK-HDF (``pharevtkhdf``)
   * - Other
     - a :term:`restart` folder ``checkpoints`` is set but no restart time, so no restart file is written; :term:`load balancing` (auto)
     - restart folder ``checkpoints`` with ``restart_time="auto"`` (resumes from it if a restart is there; none is written); :term:`load balancing` every 1000 steps

.. literalinclude:: ../../../tests/functional/harris/harris_2d.py
   :language: python
   :pyobject: config

How to run
----------

Work in an empty folder (see :doc:`/examples/index`). The test suite runs both
scripts on 10 :term:`MPI <MPI>` processes.

.. code-block:: bash

   mkdir -p harris && cd harris
   PYTHONPATH=$BUILD:$PHARE/pyphare:$PHARE \
       mpirun -n 10 python3 $PHARE/tests/functional/harris/harris_2d.py
   PYTHONPATH=$BUILD:$PHARE/pyphare:$PHARE \
       mpirun -n 10 python3 $PHARE/tests/functional/harris/harris_3d.py

What it plots
-------------

The 2D script draws, at each output time, :math:`\nabla\cdot B` (at three
colour scales: it should stay at round-off level), the MPI ranks owning each
patch, the proton density, :math:`B_x, B_y, B_z`, :math:`J_z` and two
components of the proton pressure tensor. The PNG files go to
``phare_outputs/harris_2d_plots/<number of ranks>/``. The HDF5 files are
deleted when the script ends; only the plots stay. The 3D script has no
post-processing: it writes the diagnostics.

Runtime
-------

2D: about 1 hour (estimated from a 10-minute partial run on 4 processes, machine shared). 3D: long run, not timed here; it was started on 4 processes and had not
reached the end of its 1000 steps after 180 s, at about 1 s per step).
