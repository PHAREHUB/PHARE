====================================
Magnetic jump relaxation (hybrid 1D)
====================================

**Goal:** follow a smoothed jump in the magnetic field in a hybrid simulation,
for the three particle interpolation orders, and compare the results.

The density is uniform (1) and there is no initial flow. The field
:math:`B_y` is 0.125 outside and 4 inside a slab between :math:`x = 0.2L` and
:math:`0.8L`, with :math:`\tanh` edges of width 1 (:math:`B_x = B_z = 0`). The
ion thermal speed is 0.1 and the electron temperature 0.12. The two
discontinuities then evolve into waves.

Setup
-----

.. list-table::
   :header-rows: 1

   * - Item
     - Value
   * - Dimension / domain
     - 1D, 2500 cells of size 0.2 (box length 500), periodic
   * - Model
     - hybrid, one proton :term:`population`; interpolation order 1, 2 then 3 (three runs, one after the other)
   * - Electrons
     - isothermal :term:`closure`, :math:`T_e = 0.12`
   * - Time step
     - 0.005, final time 30 (6000 steps per run)
   * - Resistivity
     - hyper-resistivity 0.01
   * - :term:`AMR`
     - no, a single :term:`level`
   * - Output
     - E, B, charge density, bulk velocity, every 10 steps

.. literalinclude:: ../../../tests/functional/shock/shock.py
   :language: python
   :pyobject: config

How to run
----------

Work in an empty folder (see :doc:`/examples/index`). The test suite runs this
script on 5 :term:`MPI <MPI>` processes.

.. code-block:: bash

   mkdir -p shock && cd shock
   PYTHONPATH=$BUILD:$PHARE/pyphare:$PHARE \
       mpirun -n 4 python3 $PHARE/tests/functional/shock/shock.py

What it checks and plots
------------------------

There is no numerical assertion. For each run the script draws :math:`B_y(x)`
at every output time and, if ``ffmpeg`` is installed, makes a movie
``shock_interp<order>.mp4`` (otherwise it prints that it skips the video). It
then removes the frames and saves ``shock_By.png``: :math:`B_y` at
:math:`t = 30` for the three interpolation orders on one plot.

Runtime
-------

About 45 minutes (estimated from a 10-minute partial run on 4 processes, machine shared).
