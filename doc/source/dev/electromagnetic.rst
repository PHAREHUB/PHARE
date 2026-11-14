.. _dev-electromagnetic:

======================================
Faraday, Ampère, Ohm and the electrons
======================================

.. container:: breadcrumb

   You are here: :doc:`/dev/index` › :doc:`/dev/one_time_step` › Faraday, Ampère, Ohm

**In one sentence:** in the hybrid model, three small operators turn the ion moments
and the magnetic field into a new electric field, and that electric field into a new
magnetic field; the electrons are only a fluid that feeds Ohm's law.

The physics behind these equations is on :doc:`/theory/hybridpic`. This page shows
where each equation lives in the code and what it needs from its neighbours.

Why it is built this way
------------------------

In the hybrid model the electrons have no equation of their own to advance in time.
They are a massless fluid: their momentum equation is solved for the electric field
instead (Ohm's law). So a field update is always the same short chain: advance B with
Faraday's law, get the current from Ampère's law, ask the electrons for their velocity
and pressure, then compute E with Ohm's law. Each link is one small class that works
on one patch.

The picture
-----------

Each stage of the hybrid time step (:term:`predictor 1, predictor 2, corrector <PPC>`) runs the same
four calls, on every patch of the level:

.. mermaid::

   flowchart TD
     F["Faraday<br/>B_new = B − Δt curl E"] --> A["Ampère<br/>J = curl B_new"]
     A --> U["Electrons::update<br/>electron velocity, pressure"]
     U --> O["Ohm<br/>E from Ve, B, Pe, J"]
     M["ion moments<br/>(charge density, charge flux)"] --> U
     click F href "electromagnetic.html#faraday" "Faraday"
     click A href "electromagnetic.html#ampere" "Ampère"
     click U href "electromagnetic.html#the-electrons" "Electrons"
     click O href "electromagnetic.html#ohm" "Ohm"
     click M href "moments.html" "Moments: deposition and sum borders"

What changes between stages is which B and E go in, and which moments the electrons
see. That order is on :doc:`/numerics/time_integration` (the predictor-corrector) and,
for the ion push in between, on :doc:`/dev/ppc_order`.

How it works
------------

#. All three operators work on one patch at a time and know nothing about AMR. They get
   a :doc:`GridLayout </dev/gridlayout>` that says where each quantity lives.
#. Derivatives come from one shared header,
   :src:`src/core/numerics/differential_operators.hpp`: ``curl``, ``grad``, ``div``,
   ``laplacian``. A derivative of a primal quantity is dual, and the other way round, so
   ``curl B`` lands exactly where ``J`` lives, and ``curl E`` where ``B`` lives.
#. When a term is needed at a place where its ingredients do not live (for example
   ``Vy × Bz`` on an ``Ex`` edge), the ingredient is averaged to that place with
   ``GridLayout::project<Source, Destination>``.
#. The level-wide loops (over patches) are thin wrappers: Faraday and Ampère in
   :src:`src/amr/solvers/solver_field_evolvers.hpp`, Ohm in
   :src:`src/amr/solvers/solver_hybrid_field_evolvers.hpp`. Ghost values are filled
   afterwards by the messengers (:doc:`/dev/messengers`).

Faraday
^^^^^^^

``Faraday::operator()`` (:src:`src/core/numerics/faraday/faraday.hpp#L28-L56`) computes
:math:`\mathbf{B}_{new} = \mathbf{B} - \Delta t\,\nabla\times\mathbf{E}` on the domain
nodes of each B component (``GridLayout::evalOnBox``). B ghosts are filled later by the
messenger.

Ampère
^^^^^^

``Ampere::operator()`` (:src:`src/core/numerics/ampere/ampere.hpp#L35-L63`) computes
:math:`\mathbf{J} = \nabla\times\mathbf{B}` (normalised units, so :math:`\mu_0 = 1`, see
:doc:`/theory/normalization`). The box it runs on is a template choice
(``AmpereBox``, :src:`src/core/numerics/ampere/ampere.hpp#L14-L20`). The hybrid solver
computes J on the domain grown by one layer, exactly the layer Ohm's law reads, so J
never needs its ghosts filled by a messenger (``FieldEvolverDispatchers::ampere_box``,
:src:`src/amr/solvers/solver_field_evolvers.hpp#L108-L120`). The MHD solver keeps the class
default: the ghost box minus its outermost layer.

Ohm
^^^

``Ohm::operator()`` (:src:`src/core/numerics/ohm/ohm.hpp#L52-L63`) computes, on the
domain edges of each E component,

.. math::

   \mathbf{E} = -\mathbf{V}_e\times\mathbf{B} - \frac{\nabla P_e}{n}
                + \eta\,\mathbf{J} - \nu\,\nabla^2\mathbf{J}

``ideal_``, ``pressure_``, ``resistive_`` and ``hyperresistive_`` are one function per
term. ``n`` is the ion charge density. ``∇Pe`` is computed once per patch and then divided
by ``n`` averaged to the edge (``pressure_``,
:src:`src/core/numerics/ohm/ohm.hpp#L173-L210`).

The electrons
^^^^^^^^^^^^^

``Electrons::update`` (:src:`src/core/data/electrons/electrons.hpp#L264-L272`) does two
things before Ohm runs:

* electron velocity, from :term:`quasineutrality` and the current
  (``StandardHybridElectronFluxComputer::computeBulkVelocity``,
  :src:`src/core/data/electrons/electrons.hpp#L70-L98`):

  .. math::

     \mathbf{V}_e = \frac{\sum_s q_s \mathbf{F}_s - \mathbf{J}}{N_e}

  where :math:`\mathbf{F}_s` is the particle flux of population :math:`s` and
  :math:`N_e = \sum_s q_s n_s` the ion charge density. All three live on the primal
  nodes; J is averaged to the nodes first.
* electron pressure, from the :term:`closure`. The only closure is isothermal,
  :math:`P_e = N_e T_e` with a constant :math:`T_e`
  (``IsothermalElectronPressureClosure``,
  :src:`src/core/data/electrons/electrons.hpp#L112-L178`). The user sets it with
  ``ElectronModel(closure="isothermal", Te=...)`` (:doc:`/user_guide/electrons`).

The tricks
----------

Charge flux, not ion bulk velocity
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** The textbook :math:`\mathbf{u}_e = \mathbf{u}_i - \mathbf{J}/(en)` uses the
ion bulk velocity. In PHARE the ion bulk velocity is mass-weighted, so with several
species of different charge-to-mass ratio it gives the wrong electron velocity.

**Fix.** The electron velocity is built from the charge flux of every population and the
charge density, as in the formula above
(:src:`src/core/data/electrons/electrons.hpp#L75-L76`). For one population of protons it
reduces to the textbook formula.

**History.** :issue:`1466` (electron velocity from charge).

Switching terms off without testing in the loop
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** Resistivity :math:`\eta` and hyper-resistivity :math:`\nu` are runtime
inputs and often zero. Testing them in every cell of the hot loop costs time; computing a
Laplacian that is multiplied by zero costs more.

**Fix.** ``Ohm::operator()`` turns the two runtime flags (``eta > 0``, ``nu > 0``) into
compile-time template arguments once per call, with ``Constexprifier`` (a small helper
that makes the compiler build one version of the loop per combination of flags)
(:src:`src/core/numerics/ohm/ohm.hpp#L56-L62`). The per-cell code then keeps or drops each
term with ``if constexpr``.

Hyper-resistivity that follows the mesh
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** Hyper-resistivity is there to damp grid-scale noise. With a constant
:math:`\nu`, its damping relative to the waves it should damp changes from one AMR level
to the next. An earlier spatial law, :math:`\nu|\mathbf{B}|/n`, followed the field but went
to zero where B vanishes, at reconnection X points, exactly where dissipation is needed.

**Fix.** With ``hyper_mode="spatial"``, the coefficient is
(``spatial_hyperresistive_``, :src:`src/core/numerics/ohm/ohm.hpp#L259-L289`)

.. math::

   \nu_{\text{eff}} = \nu\left(\frac{|\mathbf{B}|}{n + 0.1} + 1\right)\frac{1}{4^{\ell}}

with :math:`\ell` the level number. With a refinement ratio of 2, :math:`4^{-\ell}`
follows :math:`\Delta x^2`. The ``+ 1`` keeps a floor where B goes to zero; the ``0.1``
avoids a blow-up where the density is very low. The default ``hyper_mode="constant"``
uses :math:`\nu` everywhere.

**History.** :pr:`911` (safer hyper-resistivity).

.. note::

   The MHD solver has its own ``spatial`` law, :math:`\nu\,\Delta x^2(|\mathbf{B}|/\rho+1)`
   with :math:`\rho` the mass density and :math:`\Delta x` the smallest mesh size
   (:src:`src/core/numerics/constrained_transport/dissipative_electric_field.hpp#L84-L104`),
   see :doc:`/dev/mhd_solver`.

.. admonition:: Known issue
   :class: note

   In 2D, the isothermal pressure term :math:`-\nabla P_e/n` as written above can create a
   magnetic field out of density noise, because its discrete curl is not exactly zero.
   See :issue:`1527` (pressure term creates B); a fix is proposed in :pr:`1528`.

Go deeper
---------

* :doc:`/dev/gridlayout`: centerings, ``project``, the boxes the loops run on.
* :doc:`/dev/fields`: what ``Field`` and ``VecField`` are.

Related
-------

* :doc:`/dev/one_time_step` (up): where these calls sit in a time step.
* :doc:`/dev/ppc_order`: the ion push between the field stages.
* :doc:`/dev/moments`: where N_e and F_s come from.
* :doc:`/dev/messengers`: how E, B and J ghosts are filled after each stage.
* :doc:`/dev/mhd_solver`: the MHD counterpart.
