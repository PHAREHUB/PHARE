====================================
Convergence order of the MHD schemes
====================================

**Goal:** measure how fast the error of the :term:`MHD` solver drops when the
grid is refined, for each spatial reconstruction, and compare with the
expected order.

A small-amplitude Alfvén wave (:math:`10^{-6}`) travels along a uniform field
:math:`B_x = 1` in a plasma of density 1 and pressure 0.1. The Alfvén speed is 1
and the box length is 1, so after a final time of 1 the wave is back at its
starting point: the error is the mean absolute difference between
:math:`B_y` at the end and at the start. The slope of the error against the
cell size in a log-log plot is the measured order. The time stepper is ``SSPRK4_5``, chosen in the script so that the time error
does not hide the spatial one.

Two scripts do this: ``mhd_convergence`` (1D) and
``mhd_multidimensional_convergence`` (3D).

Setup, 1D (``convergence.py``)
------------------------------

.. list-table::
   :header-rows: 1

   * - Item
     - Value
   * - Dimension / domain
     - 1D, box length 1, periodic; 50, 100, 200, 400, 800 cells
   * - Model
     - ideal MHD, :math:`\gamma = 5/3`
   * - Reconstructions and expected order
     - ``Constant`` 1, ``Linear`` (VanLeer limiter) 2, ``WENO3`` 3, ``WENOZ`` 5, ``MP5`` 5
   * - Scheme
     - Rusanov :term:`Riemann solver`, ``SSPRK4_5``
   * - Time step
     - 5e-4, final time 1
   * - :term:`AMR`
     - no, a single :term:`level`
   * - Pass if
     - the fitted slope is within 15% of the expected order

.. literalinclude:: ../../../tests/functional/mhd_convergence/convergence.py
   :language: python
   :pyobject: config

Setup, 3D (``multidimensional_convergence.py``)
-----------------------------------------------

The 3D script uses an oblique wave (the mode is chosen in the script:
``Alfven`` by default, or ``Fast``, ``Slow``, ``Entropy``) on grids of
:math:`2N \times N \times N` cells for :math:`N = 16, 32, 64, 128` (box
:math:`3 \times 1.5 \times 1.5`), with the reconstructions ``Linear``,
``WENO3``, ``WENOZ`` and ``MP5``. The expected order is 2 for all of them,
because the 3D scheme is second order overall (as the script comments). The
time step is 0.0007. As shipped, the line that runs the simulation
(``Simulator(sim).run().reset()``) is commented out in this script: it only
analyses output already present in ``phare_outputs/convergence``.

How to run
----------

Work in an empty folder (see :doc:`/examples/index`). The test suite runs
both scripts on 4 :term:`MPI <MPI>` processes.

.. code-block:: bash

   mkdir -p mhd_convergence && cd mhd_convergence
   PYTHONPATH=$BUILD:$PHARE/pyphare:$PHARE \
       mpirun -n 4 python3 $PHARE/tests/functional/mhd_convergence/convergence.py

What it checks and plots
------------------------

For each reconstruction the script runs every resolution, fits the slope,
saves ``phare_outputs/convergence/convergence_<reconstruction>.png`` and
asserts that the slope is within 15% of the expected order.

Runtime
-------

On 4 :term:`MPI <MPI>` processes the 1D script passed ``Constant``, ``Linear``
and ``WENO3`` (89 s in total, other runs using the machine at the same time),
then stopped at ``WENOZ`` with ``AssertionError: WENOZ: got nan, expected 5.0``
on this build. The 3D script was checked in dry-run mode only.
