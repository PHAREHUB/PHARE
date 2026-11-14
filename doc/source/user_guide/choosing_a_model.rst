===============================
Choosing a model: hybrid or MHD
===============================

**Goal:** decide whether your problem needs kinetic ions (the
:term:`hybrid model`) or whether a single fluid (:term:`MHD`) is enough, and
know what changes in the input script.

PHARE runs one of two models. You choose with ``Simulation(model_options=...)``:

* **Hybrid** (the default, ``model_options=["HybridModel"]``): ions are
  macro-particles, electrons a massless fluid. Set up with
  ``MaxwellianFluidModel`` and ``ElectronModel``.
* **MHD** (``model_options=["MHDModel"]``): the whole plasma is one fluid,
  described by its density, velocity, pressure and magnetic field. Set up with
  ``MHDModel``.

.. important::

   One run uses one model: hybrid and MHD **cannot be combined in one run
   yet**. Running MHD on the coarse mesh levels and hybrid on the fine ones
   does not work, even though ``Simulation`` accepts
   both names in ``model_options`` and a ``max_mhd_level`` keyword without
   complaining.


Which one?
----------

* Ion kinetic effects matter (beams, non-Maxwellian ions, cyclotron or Landau
  damping, ion-scale structure): :doc:`hybrid <hybrid_inputs>`.
* Scales near the ion inertial length, without ion kinetics: MHD with the Hall
  term (``hall=True``, :doc:`mhd_inputs`), or hybrid.
* Larger scales only: :doc:`MHD <mhd_inputs>`.

.. list-table::
   :header-rows: 1
   :widths: 24 38 38

   * -
     - Hybrid
     - MHD
   * - Ions
     - particles: any velocity distribution, several populations
     - one fluid, local thermal equilibrium
   * - Electrons
     - massless fluid, isothermal closure (:doc:`electrons`)
     - part of the single fluid; ion-scale effects with ``hall=True``
   * - Smallest scales
     - ion inertial length and ion gyroradius, resolved by the grid
     - no built-in scale (ideal MHD); the ion inertial length with the Hall
       term
   * - Cost
     - grows with the number of particles; noise falls as
       :math:`1/\sqrt{\text{particles per cell}}`
     - no particles, no particle noise; much cheaper per cell
   * - Time step
     - limited by the fastest wave on the grid, the whistler, which forces
       a smaller time step when cells get smaller
     - limited by the fast magnetosonic speed (and by whistlers with the
       Hall term)
   * - Boundaries
     - periodic only
     - periodic or not (:doc:`boundary_conditions`)
   * - Input blocks
     - ``MaxwellianFluidModel`` + ``ElectronModel`` (:doc:`hybrid_inputs`)
     - ``MHDModel`` (:doc:`mhd_inputs`)
   * - Fluid outputs
     - ``FluidDiagnostics`` (ion moments)
     - ``MHDDiagnostics``

Both models write the electromagnetic field with ``ElectromagDiagnostics``,
use the same :doc:`normalized units </theory/normalization>` and the same
adaptive mesh refinement (:doc:`amr`).


What changes in the script
--------------------------

A hybrid script (as in :doc:`/getting_started/first_simulation`) needs no
``model_options``: hybrid is the default.

.. code-block:: python

   sim = ph.Simulation(cells=200, dl=0.5, time_step=0.01, final_time=50.0)
   ph.MaxwellianFluidModel(bx=1.0, protons={"density": 1.0, "vthx": 0.1,
                                            "vthy": 0.1, "vthz": 0.1})
   ph.ElectronModel(closure="isothermal", Te=0.01)

An MHD script sets ``model_options`` and the numerical scheme in
``Simulation``, then gives the fluid in ``MHDModel``; no ``ElectronModel``.
The complete, runnable version, with the list of schemes, is in
:doc:`mhd_inputs`.

.. code-block:: python

   sim = ph.Simulation(
       cells=200, dl=0.5, time_step=0.01, final_time=50.0,
       model_options=["MHDModel"],
       reconstruction="Linear", limiter="VanLeer", riemann="Rusanov",
       mhd_timestepper="TVDRK2",
   )
   ph.MHDModel(density=1.0, vx=0.0, vy=0.0, vz=0.0,
               bx=1.0, by=0.0, bz=0.0, p=0.1)


Common mistakes
---------------

* ``ValueError: Error: MHDModel requires a non-empty mhd_timestepper`` — an MHD
  run needs the scheme keywords in ``Simulation`` (see :doc:`mhd_inputs`).
* Putting both ``"HybridModel"`` and ``"MHDModel"`` in ``model_options`` to get
  a coupled run — not supported yet (see the box at the top).
* A hybrid run with cells much larger than the ion inertial length — the run
  works, but the ion kinetic physics the hybrid model is there for is not
  resolved. If the physics happens at larger scales, MHD is cheaper.
