=======================================
Energy conservation in a uniform plasma
=======================================

**Goal:** see how well the hybrid model conserves kinetic energy in a
uniform, quiet plasma, as a function of the ion thermal speed and of the
resolution.

The plasma has density 1, no bulk flow and a uniform field :math:`B_0 = 1`
along :math:`x`. Nothing should happen, so any drift of the total ion kinetic
energy measures a numerical error. The script scans the thermal speed
:math:`v_{th}` over 0.01, 0.05, 0.1, 0.3, 0.5, 0.75, 1 and 2, at two
resolutions.

Setup
-----

.. list-table::
   :header-rows: 1

   * - Item
     - Value
   * - Dimension / domain
     - 1D, periodic: 100 cells of size 0.2 (25 000 steps, time step 0.002) and 200 cells of size 0.1 (100 000 steps, time step 5e-4)
   * - Model
     - hybrid, one proton :term:`population` (100 particles per cell, default), interpolation order 1 (default)
   * - Electrons
     - isothermal :term:`closure`, :math:`T_e = 0.1`
   * - Final time
     - 50
   * - :term:`AMR`
     - no, a single :term:`level`
   * - Output
     - B and the particles (domain and level ghosts), every 50 steps

That makes 16 runs, one per pair (:math:`v_{th}`, resolution), each written in
its own folder ``vth<value>dx<value>``.

.. literalinclude:: ../../../tests/functional/conservation/conserv.py
   :language: python
   :pyobject: uniform

How to run
----------

Work in an empty folder (see :doc:`/examples/index`). The test suite runs this
script on 5 :term:`MPI <MPI>` processes.

.. code-block:: bash

   mkdir -p conservation && cd conservation
   PYTHONPATH=$BUILD:$PHARE/pyphare:$PHARE \
       mpirun -n 4 python3 $PHARE/tests/functional/conservation/conserv.py

What it checks and plots
------------------------

For each run the script sums the ion kinetic energy over the particles and the
magnetic energy over the grid, at each output time. It averages the kinetic
energy between :math:`t = 3` and 4 as the reference (to skip an initial
transient) and saves ``K.png``: the relative change of kinetic energy at the
end of the run, in percent, against :math:`v_{th}`, for both resolutions. There
is no assertion.

Runtime
-------

About 2 hours (estimated from a 10-minute partial run on 4 processes, machine shared).
