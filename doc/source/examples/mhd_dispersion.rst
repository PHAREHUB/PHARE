=============================
Whistler waves (Hall MHD, 1D)
=============================

**Goal:** propagate small whistler waves in 1D with the Hall :term:`MHD`
solver, at a fine and a coarse resolution, and write the field for a
dispersion analysis.

The plasma has density 1, pressure 1, a field :math:`B_0 = 1` along
:math:`x`, and four transverse circularly polarized modes (numbers 1, 2, 4, 8,
random phases, amplitude 0.01 in the field). The script defines the wave
frequency of mode ``m`` as
:math:`\omega = \tfrac{k^2}{2}\left(\sqrt{1 + 4/k^2} + 1\right)` with
:math:`k = 2\pi m/L_x`, and runs for 10 periods of that frequency.

Setup
-----

.. list-table::
   :header-rows: 1

   * - Item
     - ``high`` resolution
     - ``low`` resolution
   * - Cells / size
     - 128 / 0.05
     - 128 / 0.8
   * - Time step
     - 0.0006
     - 0.077
   * - Model
     - Hall MHD (``hall=True``), :math:`\gamma = 5/3`
     - same
   * - Scheme
     - linear reconstruction, VanLeer limiter, Rusanov :term:`Riemann solver`, ``TVDRK2``
     - same
   * - :term:`AMR`
     - no, a single :term:`level`
     - no
   * - Output
     - B at every time step
     - B at every time step

.. literalinclude:: ../../../tests/functional/mhd_dispersion/whistler1d.py
   :language: python
   :pyobject: config

How to run
----------

Work in an empty folder (see :doc:`/examples/index`). The test suite runs this
script on 4 :term:`MPI <MPI>` processes. It runs the ``high`` case, then the
``low`` case.

.. code-block:: bash

   mkdir -p mhd_dispersion && cd mhd_dispersion
   PYTHONPATH=$BUILD:$PHARE/pyphare:$PHARE \
       mpirun -n 4 python3 $PHARE/tests/functional/mhd_dispersion/whistler1d.py

What it checks
--------------

Nothing: the script only runs the simulations and leaves the B files in
``phare_outputs/whistler/high`` and ``phare_outputs/whistler/low``. Measuring
the frequencies is left to you (the hybrid counterpart,
:doc:`dispersion`, does it).

Runtime
-------

Long run: not timed here (cluster recommended). It had not finished the
``high`` case after 400 s on 4 processes.
