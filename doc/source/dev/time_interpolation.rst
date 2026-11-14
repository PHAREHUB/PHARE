.. _dev-time-interpolation:

===========================================
Interpolation in time: ghosts between steps
===========================================

.. container:: breadcrumb

   You are here: :doc:`Developer guide </dev/index>` › :doc:`Sub-cycling </dev/subcycling>` › Interpolation in time

**Applies to:** hybrid ✓ | MHD ✓ (fields only; the particle side is on :doc:`level_ghost_particles`)

**In one sentence:** a fine level-ghost cell gets, at the time of a fine
sub-step, a straight-line mix of two coarse values: the one saved at the start of
the coarse step and the one already computed for its end.

Why
---

The coarse level steps first (:doc:`subcycling`). When the fine level needs its
level-ghost cells, the coarse level is already at :math:`t+\Delta t`, but the fine
sub-step lives at some time in between. The coarse level therefore keeps a copy
of its values at :math:`t`, and the fine level mixes the two.

.. plot::
   :include-source: false
   :alt: One coarse step, four fine sub-steps, the two saved coarse times, and the times at which B and the averaged E are filled in sub-step 2.

   from phare_draw import *
   fig, (a,) = strip(1, aspect=False, height=3.0)
   amr_clock(a, side_labels=True, shade=[2],
             coarse_dots=[(0, "saved at t"), (1, "already at t+dt")],
             marks=[(0.375, "Eavg"), (0.5, "B")])
   finish(fig, [a], align=False)

Look at sub-step 2 (shaded): the averaged electric field :math:`E_\text{avg}` is
filled in its middle, :math:`t + 3\Delta t/8`, and B at its end,
:math:`t + \Delta t/2`. Each field has its own fill time, so each has its own
weight.

The two coarse values
---------------------

* **Old** is a copy the coarse level makes at the start of each of its own steps
  (``prepareStep``), stamped with that time. Hybrid saves B and E; MHD saves
  :math:`\rho`, :math:`\rho \mathbf v`, :math:`E_\text{tot}` and B
  (:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L580-L597`,
  :src:`src/amr/messengers/mhd_messenger.hpp#L324-L342`).

.. claim C-ti-01

* **New** is the coarse model array itself, already advanced to the end of the
  coarse step when the fine level steps
  (:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L789-L831`,
  :src:`src/amr/messengers/mhd_messenger.hpp#L554-L590`).

.. claim C-ti-02

* The mix is :math:`(1-\alpha)\,\text{old} + \alpha\,\text{new}` with
  :math:`\alpha = (t_\text{fill} - t_\text{old})/(t_\text{new} - t_\text{old})`
  (``FieldLinearTimeInterpolate``,
  :src:`src/amr/data/field/time_interpolate/field_linear_time_interpolate.hpp#L38`,
  :src:`src/amr/data/field/time_interpolate/field_linear_time_interpolate.hpp#L95-L98`).

.. claim C-ti-03

The weight, sub-step by sub-step
--------------------------------

With four equal sub-steps (hybrid, see :doc:`subcycling`), old at :math:`t` and
new at :math:`t+\Delta t`:

.. list-table::
   :header-rows: 1
   :widths: 14 28 28

   * - Sub-step
     - :math:`\alpha` for B, Bpred, E (filled at the sub-step end)
     - :math:`\alpha` for :math:`E_\text{avg}` (filled at the middle)
   * - 1
     - 0.25
     - 0.125
   * - 2
     - 0.5
     - 0.375
   * - 3
     - 0.75
     - 0.625
   * - 4
     - 1
     - 0.875

Always say which field a weight belongs to: 0.375 is the :math:`E_\text{avg}`
weight of sub-step 2, never the B weight.

.. claim C-ti-05

At :math:`\alpha = 1` (end of the last sub-step) only the new value is used.
SAMRAI expects :math:`t_\text{old} < t_\text{fill} \le t_\text{new}`, checked in
debug builds only (``RefineTimeTransaction::timeInterpolate``, SAMRAI
``xfer/RefineTimeTransaction.cpp``).

.. claim C-ti-11

Which fields, and when they are filled
--------------------------------------

.. list-table::
   :header-rows: 1
   :widths: 10 14 26 24 26

   * - Model
     - Field
     - Old copy
     - New
     - Filled at
   * - hybrid
     - B, Bpred
     - ``Bold_`` (messenger)
     - model B
     - end of the sub-step
   * - hybrid
     - E
     - ``Eold_`` (messenger)
     - model E
     - end of the sub-step
   * - hybrid
     - :math:`E_\text{avg}`
     - ``Eold_``
     - model E (not the coarse :math:`E_\text{avg}`)
     - middle of the sub-step
   * - MHD
     - :math:`\rho`, :math:`\rho \mathbf v`, :math:`E_\text{tot}`, B
     - ``rhoOld_``, ``rhoVold_``, ``EtotOld_``, ``Bold_``
     - model arrays
     - TVDRK2: every fill of the step is at the step end; other integrators not described
   * - MHD
     - E
     - not taken from the coarse level
     - —
     - —

.. claim C-ti-06 C-ti-07 C-ti-08 C-ti-10

* Hybrid fill times are in the PPC step: B and Bpred at ``newTime``, E at
  ``newTime``, :math:`E_\text{avg}` at ``(currentTime + newTime)/2``
  (:src:`src/amr/solvers/solver_ppc.hpp#L427`,
  :src:`src/amr/solvers/solver_ppc.hpp#L508`,
  :src:`src/amr/solvers/solver_ppc.hpp#L528`,
  :src:`src/amr/solvers/solver_ppc.hpp#L561`).

.. claim C-ti-04

* The hybrid current J is never communicated: each level computes it from its
  time-interpolated B (:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L834-L835`).

.. claim C-ti-09

Time first, space second
------------------------

.. flow::

   two coarse values | old (saved at t) and new (computed, at t+dt)
   mixed on the coarse grid
   coarse values at the fill time | one straight-line mix
   then refined in space
   fine ghost value | level-ghost cell | /dev/field_refinement

SAMRAI builds the coarse values at the fill time first, in a recursive coarse
schedule that uses time interpolation, and then refines them in space
(``RefineSchedule::recursiveFill``, SAMRAI ``xfer/RefineSchedule.cpp``).

.. claim C-ti-12

Remark, for fields refined with the composite rule (all but B): both steps are
linear in the coarse values, so the order would not change the result. Not
claimed for B, whose second stage also keeps faces already written on the fine
level (see :doc:`field_refinement`).

.. claim C-fr-16

Example, 1D
-----------

On the running example (:doc:`field_refinement`) the level-ghost cell 3 of
:math:`B_y` is 1.275 when refined from the coarse values at :math:`t`, and 1.32
when refined from those at :math:`t+\Delta t`. In sub-step 2, B is filled at
:math:`\alpha = 0.5`: :math:`0.5 \times 1.275 + 0.5 \times 1.32 = 1.2975`.

.. claim C-ex-04 C-ti-05

Regridding and a new level
--------------------------

The first coarse step of a new level, and what a regrid does to these copies,
are not described here: see :doc:`tagging_regridding`.

Next: :doc:`refluxing`, how the coarse level takes the fine answer when the levels meet.

Go deeper
---------

* :doc:`subcycling`: why the coarse level runs ahead.
* :doc:`field_refinement`: the space part.
* :doc:`level_ghost_particles`: the same idea for particles.
* :doc:`messengers`: where these refiners are registered.

Related
-------

* :doc:`ghost_cells`: patch ghost cells versus level ghost cells.
* :doc:`/numerics/amr`: the AMR algorithm for users.
