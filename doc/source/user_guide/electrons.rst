=========
Electrons
=========

**Goal:** set the electron fluid of a hybrid run, and the dissipation terms of
Ohm's law that go with it.

In the :term:`hybrid model`, ions are particles but electrons are a massless
fluid. PHARE never moves electrons: it uses their momentum equation to compute
the electric field. What you choose is the electron :term:`closure` (how the
electron pressure is computed) and, in the ``Simulation`` block, the
resistivity and hyper-resistivity. An MHD run does not use this page.


Minimal example
---------------

.. code-block:: python

   import pyphare.pharein as ph

   sim = ph.Simulation(
       # ... grid and time, as in the first simulation ...
       resistivity=0.001,          # eta, optional
       hyper_resistivity=0.01,     # nu, optional
       hyper_mode="spatial",       # optional
   )

   # ... ph.MaxwellianFluidModel(...) ...

   ph.ElectronModel(closure="isothermal", Te=0.01)

``ElectronModel`` is required in every hybrid run. The three ``Simulation``
keywords are optional (defaults below).


What PHARE computes from the electrons
--------------------------------------

At every step, PHARE builds the electron fluid from the ions and the magnetic
field, then gets :math:`\mathbf{E}` from Ohm's law (in
:doc:`normalized units </theory/normalization>`):

.. mermaid::

   flowchart TD
     I["ion particles"] --> M["ion moments<br/>densities and fluxes"]
     M --> N["electron density<br/>n_e = sum of q_s n_s"]
     B["magnetic field B"] --> J["current J = curl B"]
     M --> V["electron velocity<br/>V_e = (sum of q_s n_s V_s - J) / n_e"]
     J --> V
     N --> P["electron pressure<br/>P_e = n_e T_e"]
     V --> E["Ohm's law gives E"]
     P --> E
     J --> E

* **Density:** the plasma is quasineutral, so the electron density is the ion
  charge density, :math:`n_e = \sum_s q_s n_s` over all ion populations.
* **Velocity:** electrons carry the part of the current the ions do not,
  :math:`\mathbf{V}_e = \left(\sum_s q_s n_s \mathbf{V}_s - \mathbf{J}\right)/n_e`,
  with :math:`\mathbf{J} = \nabla\times\mathbf{B}`.
* **Pressure:** with the isothermal closure,
  :math:`P_e = n_e T_e`, with :math:`T_e` the same everywhere and at all
  times.
* **Ohm's law:**

  .. math::

     \mathbf{E} = -\mathbf{V}_e\times\mathbf{B} - \frac{\nabla P_e}{n_e}
                  + \eta\,\mathbf{J} - \nu\,\nabla^2\mathbf{J}

  The first two terms are the ideal and pressure terms. The last two are the
  resistive term (:math:`\eta`, ``resistivity``) and the hyper-resistive term
  (:math:`\nu`, ``hyper_resistivity``).

The equations and how they are solved are in :doc:`/theory/hybridpic`.


Why hyper-resistivity
---------------------

In a hybrid code, the shortest waves on the grid (whistler waves) are the
fastest, and they collect numerical noise at the scale of one cell. The
hyper-resistive term :math:`-\nu\nabla^2\mathbf{J}` damps the smallest scales
much more than the large ones, so it removes this noise while leaving the
physics at larger scales almost untouched. This is why it is on by default,
with a small value.

``hyper_mode`` chooses how :math:`\nu` varies:

* ``"constant"``: :math:`\nu` is ``hyper_resistivity`` everywhere.
* ``"spatial"``: the coefficient depends on the local field and density, and
  on the mesh :term:`level` :math:`\ell` (0 for the coarsest):

  .. math::

     \nu_{\text{eff}} = \nu\,\left(\frac{|\mathbf{B}|}{n_e + 0.1} + 1\right)\frac{1}{4^{\ell}}

  It damps more where the whistler speed, which grows with :math:`|\mathbf{B}|/n`,
  is large, and less on finer levels, whose cells are smaller (the factor
  :math:`1/4^\ell` follows the square of the cell size, as :math:`\nabla^2` does).


Options
-------

``ElectronModel``:

.. list-table::
   :header-rows: 1
   :widths: 25 55 20

   * - Option
     - Meaning
     - Default
   * - ``closure``
     - how the electron pressure is computed; only ``"isothermal"`` exists
     - required
   * - ``Te``
     - electron temperature, uniform and constant (pressure = ``Te`` times
       density)
     - 0.1

In ``Simulation``:

.. list-table::
   :header-rows: 1
   :widths: 25 55 20

   * - Option
     - Meaning
     - Default
   * - ``resistivity``
     - :math:`\eta`, resistive term :math:`\eta\mathbf{J}`
     - 0
   * - ``hyper_resistivity``
     - :math:`\nu`, hyper-resistive term :math:`-\nu\nabla^2\mathbf{J}`
     - 1e-4
   * - ``hyper_mode``
     - ``"constant"`` or ``"spatial"`` (see above)
     - ``"constant"``

``Te`` is in the units of :doc:`/theory/normalization` (temperatures in
:math:`m_p v_0^2`), like the ion temperature :math:`m v_{th}^2`: with
``Te=0.01`` and ion thermal velocity 0.1, protons and electrons have the same
temperature.

MHD runs have their own dissipation keywords (``eta``, ``nu``): see
:doc:`mhd_inputs`.


Common mistakes
---------------

* ``RuntimeError: Error - no electrons registered to this Simulation`` — a
  hybrid run needs ``ph.ElectronModel(...)``.
* ``AttributeError: 'NoneType' object has no attribute 'dict_path'`` (inside an
  error from ``Simulator.setup()``) — ``closure`` is not ``"isothermal"``, the
  only closure that exists. Check the spelling.
* ``ValueError: Error: hyper_mode should be 'constant' or 'spatial'`` — the
  value is case-sensitive.
* Grid-scale wiggles in E or B that grow with time — the hyper-resistivity may
  be too small for your cell size and time step. Increase ``hyper_resistivity``
  a little, or try ``hyper_mode="spatial"``.


Reference
---------

.. autoclass:: pyphare.pharein.ElectronModel
