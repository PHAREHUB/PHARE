=================================
Ion-ion beam instability (hybrid)
=================================

**Goal:** measure the growth rate of the right-hand ion-ion beam instability
and compare it to the expected value, 0.09.

A dense background of protons (density 1, at rest) is crossed by a tenuous
proton beam (density 0.01, velocity 5 along the field), both with thermal speed
:math:`\sqrt{0.1}`, in a uniform field :math:`B_0 = 1` along :math:`x`. The
box length, 33, is one wavelength of the most unstable mode
(:math:`k = 0.19`, as noted in the script), so the instability grows in
Fourier mode 1 of the transverse field.

Setup
-----

.. list-table::
   :header-rows: 1

   * - Item
     - Value
   * - Dimension / domain
     - 1D, 165 cells of size 0.2 (box length 33), periodic
   * - Model
     - hybrid, two proton :term:`populations <population>`: ``main`` and ``beam``
   * - Electrons
     - isothermal :term:`closure`, :math:`T_e = 0`
   * - Time step
     - 0.001, final time 100 (100 000 steps)
   * - Resistivity
     - hyper-resistivity 0.01
   * - :term:`AMR`
     - yes: :term:`tagging`, up to 3 :term:`levels <level>`, :term:`patch` sizes 20 to 60
   * - Output
     - B, E and the particles of both populations, every 0.1

.. literalinclude:: ../../../tests/functional/ionIonBeam/ion_ion_beam1d.py
   :language: python
   :pyobject: config

How to run
----------

Work in an empty folder (see :doc:`/examples/index` for ``$PHARE`` and
``$BUILD``). The test suite runs this script on 2 :term:`MPI <MPI>` processes.

.. code-block:: bash

   mkdir -p ion_ion_beam && cd ion_ion_beam
   PYTHONPATH=$BUILD:$PHARE/pyphare:$PHARE \
       mpirun -n 2 python3 $PHARE/tests/functional/ionIonBeam/ion_ion_beam1d.py

What it checks and plots
------------------------

The script takes the Fourier mode 1 of :math:`B_y - iB_z` at each output time
and fits an exponential :math:`a\,e^{\gamma t}` on the growing phase. It
asserts :math:`|\gamma - 0.09| < 0.02`. The figure ``ion_ion_beam1d.png`` shows
the mode amplitude with the fit, and, at three times, the ion distribution in
:math:`(x, v_x)` with :math:`E_y(x)`.

Runtime
-------

About 20 minutes (estimated from a 10-minute partial run on 4 processes, machine shared).
