=========================================================
Advected tangential discontinuities, with and without AMR
=========================================================

**Goal:** advect two 1D tangential discontinuities with a flowing plasma, once
on a single grid and once with :term:`AMR`, and check that refinement does not
break the plasma state.

A tangential discontinuity separates two plasmas in total-pressure balance. Here
:math:`B_y` goes from -1 to 1 at :math:`x = L/4` and back to -1 at
:math:`x = 3L/4`, with :math:`B_z = 0.5` and :math:`B_x = 0`. The density is
uniform (1) and the temperature follows from pressure balance,
:math:`T = (1 - B^2/2)/n`. The whole plasma flows at :math:`v_x = 2`, so the
profile is carried along without deformation in the exact solution.

Setup
-----

.. list-table::
   :header-rows: 1

   * - Item
     - Value
   * - Dimension / domain
     - 1D, 200 cells of size 1 (box length 200), periodic
   * - Model
     - hybrid, one proton :term:`population` (100 particles per cell), interpolation order 1
   * - Electrons
     - isothermal :term:`closure`, :math:`T_e = 0.12`
   * - Time step
     - 0.04, final time 11 (one output, at the end)
   * - Resistivity
     - hyper-resistivity 0.01
   * - :term:`AMR`
     - two runs: ``noRefinement`` (single :term:`level`) and ``withTagging`` (:term:`tagging`, up to 3 levels)
   * - Output
     - E, B, charge density, bulk velocity and protons, at the final time

.. literalinclude:: ../../../tests/functional/tdtagged/td1dtagged.py
   :language: python
   :pyobject: config

How to run
----------

Work in an empty folder (see :doc:`/examples/index`). The test suite runs this
script on 4 :term:`MPI <MPI>` processes (with 4 processes, level 0 has patch
boundaries that the discontinuities cross, which exercises tagging across
:term:`patches <patch>`).

.. code-block:: bash

   mkdir -p tdtagged && cd tdtagged
   PYTHONPATH=$BUILD:$PHARE/pyphare:$PHARE \
       mpirun -n 4 python3 $PHARE/tests/functional/tdtagged/td1dtagged.py

What it checks and plots
------------------------

It saves ``tdtagged1d.png``: the patches of levels 0 to 2, and :math:`B_y` and :math:`J_z` for both runs, compared
to the initial profile moved by :math:`V t`, with zooms on the two
discontinuities. It then raises an error if one of these checks fails:

* the right-hand sheet of the tagged run, fitted with a tanh, has thickness
  1 ± 0.5 and sits at :math:`3L/4 + V t` ± 0.5;
* every particle of the refined levels lies inside its :term:`patch` (the
  domain-particle check of PHARE's advance tests);
* on level 0 the number of particles equals ``cells * 100``;
* on level 0 the mean :math:`B_y` is 0 and the mean :math:`B_z` is 0.5, to
  round-off.

Runtime
-------

77 s on 4 :term:`MPI <MPI>` processes (other runs were using the same machine
at the same time).
