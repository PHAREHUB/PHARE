==============================================================
Hybrid dispersion relation (whistler and Alfvén/ion-cyclotron)
==============================================================

**Goal:** check that the hybrid model propagates small waves along the field at
the frequency given by the dispersion relation, on the whistler (right-hand)
and the Alfvén/ion-cyclotron (left-hand) branches.

The plasma is cold: ions with a very small thermal speed (:math:`10^{-3}`),
massless electrons with :math:`T_e = 0`, and a uniform field :math:`B_0 = 1`
along :math:`x`. The script loads small-amplitude (:math:`10^{-2}`) waves of
both branches at once, with random phases, then measures their frequencies
:math:`\omega` from the output and compares them with the theory. Units are the
ion inertial length and the inverse ion cyclotron frequency; the dispersion
relation is :math:`\omega^2 + k^2 d_i\,\omega - k^2 = 0` for the two branches.
The full derivation, figures and the discussion of the scheme's own dispersion
are in :doc:`/validation/hybrid_dispersion`.

Setup
-----

.. list-table::
   :header-rows: 1

   * - Item
     - Value
   * - Dimension / domain
     - 1D, periodic. ``low_k``: 128 cells of size 0.8, modes 1, 2, 4, 8, time step 0.02, final time 100. ``high_k``: 64 cells of size 0.1, modes 1 to 4, time step 0.001, final time 10
   * - Model
     - hybrid, one proton :term:`population`; 100 particles per cell (``low_k``), 400 (``high_k``)
   * - :term:`Interpolation order <interpolation order>`
     - 1, 2 and 3: pass it on the command line, or run all three by default
   * - Electrons
     - isothermal :term:`closure`, :math:`T_e = 0`
   * - Resistivity
     - none
   * - :term:`AMR`
     - no, a single :term:`level`
   * - Output
     - B only

.. literalinclude:: ../../../tests/functional/dispersion/dispersion.py
   :language: python
   :pyobject: config

How to run
----------

Work in an empty folder (see :doc:`/examples/index`). The optional argument is
the interpolation order (1, 2 or 3), and may be followed by a run name
(``low_k`` or ``high_k``); without it all runs are done at all orders. The test
suite runs one test per order, each on 5 :term:`MPI <MPI>` processes.

.. code-block:: bash

   mkdir -p dispersion && cd dispersion
   PYTHONPATH=$BUILD:$PHARE/pyphare:$PHARE \
       mpirun -n 4 python3 $PHARE/tests/functional/dispersion/dispersion.py 1

What it checks and plots
------------------------

For each wave it measures :math:`\omega`, prints a table (theory, grid-corrected
theory, measured, relative errors, amplitude ratio) and saves
``dispersion_o<order>.png``. It raises an error if one measured frequency is
off the grid-corrected value by more than 1%, or if a wave amplitude ends
outside 0.3 to 3 times its initial value.

Runtime
-------

Interpolation order 1 (both runs): 117 s on 4 :term:`MPI <MPI>` processes,
with other runs using the machine at the same time. The check passed.
