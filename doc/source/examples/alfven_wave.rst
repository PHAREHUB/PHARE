====================
Alfvén wave (hybrid)
====================

**Goal:** check that a circularly polarized Alfvén wave travels at the Alfvén
speed in the :term:`hybrid PIC` model, over a long run.

The plasma is uniform (density 1) in a uniform magnetic field :math:`B_0 = 1`
along :math:`x`. A small perturbation (amplitude 0.01) is added to
:math:`B_y, B_z` and to the ion velocity, with one wavelength in the box:
:math:`B_y = 0.01\cos(2\pi x/L)`, :math:`B_z = 0.01\sin(2\pi x/L)`, and the same
for :math:`v_y, v_z`. In PHARE's units the Alfvén speed is 1, so the wave must
move at phase speed :math:`V_\phi = 1`. The script measures this speed from the
simulation output.

Setup
-----

.. list-table::
   :header-rows: 1

   * - Item
     - Value
   * - Dimension / domain
     - 1D, 1000 cells of size 1 (box length 1000), periodic
   * - Model
     - hybrid, one proton :term:`population` (100 particles per cell, default; interpolation order 1, default)
   * - Electrons
     - isothermal :term:`closure`, :math:`T_e = 0`
   * - Time step
     - 0.01 (final time 1000, so 100 000 steps)
   * - Resistivity
     - hyper-resistivity 0.001
   * - :term:`AMR`
     - no, a single :term:`level`
   * - Output
     - E, B, charge density and bulk velocity, 11 times from 0 to 1000

.. literalinclude:: ../../../tests/functional/alfven_wave/alfven_wave1d.py
   :language: python
   :pyobject: config

How to run
----------

``$PHARE`` is your clone of PHARE and ``$BUILD`` your build folder
(see :doc:`/examples/index`). Work in an empty folder: the output is written
there. The test suite runs this script on 2 :term:`MPI <MPI>` processes.

.. code-block:: bash

   mkdir -p alfven_wave && cd alfven_wave
   PYTHONPATH=$BUILD:$PHARE/pyphare:$PHARE \
       mpirun -n 2 python3 $PHARE/tests/functional/alfven_wave/alfven_wave1d.py

What it checks and plots
------------------------

At each output time the script fits a cosine to :math:`B_y(x)`, which gives the
phase :math:`\phi(t)` and the wavenumber :math:`k`; the phase speed is
:math:`V_\phi = (d\phi/dt)/k`. It saves ``alfven_wave.png``: :math:`B_y` at
:math:`t = 0`, 500 and 1000, a dashed theory curve, and the mean
:math:`V_\phi` in the title. A final assertion compares the measured speed to
1 (5 % tolerance); as written, it passes as soon as one sample of
:math:`V_\phi` is within 5 % of 1, so look at the plot as well.

Runtime
-------

About 30 minutes (estimated from a 10-minute partial run on 4 processes, machine shared).

.. seealso:: :doc:`/getting_started/first_simulation` is a short version of the same problem.
