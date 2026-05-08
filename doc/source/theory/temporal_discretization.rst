
=======================
Temporal discretization
=======================

PHARE advances fields and particles in time using a second-order
predictor-corrector scheme combined with a Boris pusher for the particle
equations of motion. This page describes each algorithmic step, the stability
constraint that governs the choice of time step, and how AMR subcycling is
organised across refinement levels.

See also: :doc:`hybridpic`, :doc:`spatial_discretization`, :doc:`amr`.


The predictor-corrector scheme
-------------------------------

The hybrid equations (see :doc:`hybridpic`) couple fields and particles
non-linearly: the electric field depends on ion moments, which in turn depend
on where the particles are after being pushed with that electric field. A
single explicit Euler step would be only first-order accurate and would
introduce energy errors that grow quickly. PHARE instead uses a two-pass
predictor-corrector that achieves second-order accuracy in time.

At each step the ghost cells of every patch must be filled from neighbouring
patches or from coarser-level data via the messenger system (see :doc:`amr`)
before any field solve or push is performed. The sequence below applies to a
single AMR level; subcycling across levels is described in
`AMR subcycling`_.

**Step 1 — Predictor 1**

Advance the magnetic field from time :math:`t^n` to a predicted state
:math:`\mathbf{B}^*` using Faraday's law:

.. math::

   \mathbf{B}^* = \mathbf{B}^n - \Delta t\,\nabla \times \mathbf{E}^n

Compute a predicted current density :math:`\mathbf{J}^*` via Ampere's law
applied to :math:`\mathbf{B}^*`:

.. math::

   \mu_0\,\mathbf{J}^* = \nabla \times \mathbf{B}^*

Evaluate the generalised Ohm's law to obtain a predicted electric field
:math:`\mathbf{E}^*`:

.. math::

   \mathbf{E}^* = \mathrm{Ohm}\!\left(n^n,\,\mathbf{V}_e^n,\,P_e^n,\,
                                       \mathbf{B}^*,\,\mathbf{J}^*\right)

Ghost cells are filled from the messenger before each operator is applied.

**Step 2 — Time-centred averages**

Form averages of the predicted and current-time fields:

.. math::

   \bar{\mathbf{B}} = \frac{\mathbf{B}^n + \mathbf{B}^*}{2}, \qquad
   \bar{\mathbf{E}} = \frac{\mathbf{E}^n + \mathbf{E}^*}{2}

**Step 3 — First particle push (domain particles only)**

Push domain particles with :math:`\bar{\mathbf{B}}` and
:math:`\bar{\mathbf{E}}` using the Boris algorithm (see `The Boris pusher`_).
Deposit ion moments — number density :math:`n` and bulk flux
:math:`n\mathbf{V}_i` — onto the grid. Ghost particles are not pushed at this
stage; only the interior particle population contributes to the moment deposit.

**Step 4 — Predictor 2**

Repeat the predictor sequence with the updated ion moments:

.. math::

   \mathbf{B}^* &= \mathbf{B}^n - \Delta t\,\nabla \times \bar{\mathbf{E}} \\[4pt]
   \mu_0\,\mathbf{J}^* &= \nabla \times \mathbf{B}^* \\[4pt]
   \mathbf{E}^* &= \mathrm{Ohm}\!\left(n^{n+1/2},\,\mathbf{V}_e^{n+1/2},\,
                                        P_e^{n+1/2},\,\mathbf{B}^*,\,\mathbf{J}^*\right)

Ghost cells are refilled from the messenger before each evaluation.

**Step 5 — Update time-centred averages**

Recompute the averages using the new predictions:

.. math::

   \bar{\mathbf{B}} = \frac{\mathbf{B}^n + \mathbf{B}^*}{2}, \qquad
   \bar{\mathbf{E}} = \frac{\mathbf{E}^n + \mathbf{E}^*}{2}

**Step 6 — Second particle push (all particles)**

Push all particles — both domain and ghost — with the updated
:math:`\bar{\mathbf{B}}` and :math:`\bar{\mathbf{E}}`. The inclusion of ghost
particles ensures that particles which cross patch boundaries during the step
carry consistent velocities. Moments are deposited again to provide the
:math:`t^{n+1}` ion state.

**Step 7 — Corrector**

Advance the fields to :math:`t^{n+1}` using the fully updated moments:

.. math::

   \mathbf{B}^{n+1} &= \mathbf{B}^n - \Delta t\,\nabla \times \bar{\mathbf{E}} \\[4pt]
   \mu_0\,\mathbf{J}^{n+1} &= \nabla \times \mathbf{B}^{n+1} \\[4pt]
   \mathbf{E}^{n+1} &= \mathrm{Ohm}\!\left(n^{n+1},\,\mathbf{V}_e^{n+1},\,
                                            P_e^{n+1},\,\mathbf{B}^{n+1},\,
                                            \mathbf{J}^{n+1}\right)

Ghost cells are filled from the messenger one final time before the corrector
Ohm solve so that patch boundaries are consistent at :math:`t^{n+1}`.


The Boris pusher
----------------

Each particle is advanced individually under the Lorentz force:

.. math::

   m\,\frac{d\mathbf{v}}{dt} = q\!\left(\mathbf{E} + \mathbf{v}\times\mathbf{B}\right),
   \qquad
   \frac{d\mathbf{x}}{dt} = \mathbf{v}

The Boris algorithm splits the velocity update into two electric half-kicks
separated by a magnetic rotation. This splitting makes the rotation exactly
volume-preserving in velocity space, which prevents spurious secular energy
growth.

**Step 1 — First electric half-kick**

.. math::

   \mathbf{v}^- = \mathbf{v}^n + \frac{q}{m}\,\mathbf{E}\,\frac{\Delta t}{2}

**Step 2 — Magnetic rotation**

Define the rotation vectors:

.. math::

   \mathbf{t} = \frac{q}{m}\,\mathbf{B}\,\frac{\Delta t}{2},
   \qquad
   \mathbf{s} = \frac{2\,\mathbf{t}}{1 + |\mathbf{t}|^2}

Perform the rotation:

.. math::

   \mathbf{v}' &= \mathbf{v}^- + \mathbf{v}^- \times \mathbf{t} \\[4pt]
   \mathbf{v}^+ &= \mathbf{v}^- + \mathbf{v}' \times \mathbf{s}

**Step 3 — Second electric half-kick**

.. math::

   \mathbf{v}^{n+1} = \mathbf{v}^+ + \frac{q}{m}\,\mathbf{E}\,\frac{\Delta t}{2}

**Step 4 — Position update**

.. math::

   \mathbf{x}^{n+1} = \mathbf{x}^n + \mathbf{v}^{n+1}\,\Delta t

The electric and magnetic fields are interpolated from the Yee grid to the
particle position before the push (see :doc:`spatial_discretization`).


Stability constraint
--------------------

The explicit time integration is stable only if particles do not travel more
than one grid cell per time step. If a particle displacement exceeds
:math:`\Delta x` during a push, PHARE raises an error and halts the
simulation.

Users must choose the time step small enough to satisfy:

.. math::

   \Delta t < \frac{\Delta x}{\max|\mathbf{v}|}

where the maximum is taken over all particles and all velocity components.
In practice a safety factor of a few tenths is advisable because the maximum
particle speed can increase during the run.


AMR subcycling
--------------

In an AMR hierarchy, finer levels resolve smaller spatial scales and must
therefore advance with smaller time steps. A naive choice would follow the
CFL-like scaling :math:`\Delta t \propto \Delta x`, giving one fine substep
per coarse step for a refinement ratio of 2. However, the hybrid model
supports whistler waves, whose dispersion relation imposes the much stricter
scaling:

.. math::

   \Delta t \propto \Delta x^2

Consequently, the time step ratio between adjacent levels is the **square** of
the spatial refinement ratio:

.. math::

   \Delta t_{\mathrm{fine}} = \frac{\Delta t_{\mathrm{coarse}}}{r^2}

where :math:`r` is the refinement ratio. With the default ratio
:math:`r = 2`, each fine level performs **4 substeps** for every single coarse
step. This is more expensive than ratio-1 subcycling but is necessary to keep
the whistler branch stable on fine patches.

Between substeps, fine-level ghost cells are filled by time-interpolation from
coarse-level field data stored at the beginning and end of the coarse step.
After all substeps are complete, fine-level data are synchronised back to the
coarse level via restriction and flux correction (see :doc:`amr`).
