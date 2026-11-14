.. _dev-ppc-order:

=====================================
Push and deposition order in the PPC
=====================================

.. container:: breadcrumb

   You are here: :doc:`/dev/index` › :doc:`/dev/one_time_step` › Push and
   deposition order

**In one sentence:** the hybrid solver ``SolverPPC`` advances one :term:`level` by one
step with a :term:`predictor-corrector scheme <PPC>` in which particles are pushed twice: a
first push only *predicts* the ion moments at the new time, the second push
really moves the particles.

Why two pushes
--------------

In the hybrid model, the electric field comes from :term:`Ohm's law`, which needs the
ion density and bulk velocity. To advance the fields from time :math:`n` to
:math:`n+1` consistently, the solver needs the ion moments at :math:`n+1`, which
need the particles pushed with fields at :math:`n+1/2`, which need :math:`E` at
:math:`n+1`... The scheme breaks this loop by first *predicting* the fields and
moments, then correcting them. The physics of the scheme is on
:doc:`/numerics/time_integration`; this page is about what the code does with
the particles at each stage.

The picture
-----------

.. mermaid::

   flowchart TD
     P1["<b>predictor 1</b><br/>Faraday, Ampère, Ohm<br/>moments at n → predicted B, E"] --> A1["<b>average</b><br/>(fields at n + predicted) / 2"]
     A1 --> M1["<b>push 1</b>, moments only<br/>a copy of the particles is pushed<br/>→ predicted moments at n+1"]
     M1 --> P2["<b>predictor 2</b><br/>Faraday, Ampère, Ohm<br/>predicted moments → predicted B, E"]
     P2 --> A2["<b>average</b>"]
     A2 --> M2["<b>push 2</b>, for real<br/>particles move, change patch<br/>→ moments at n+1"]
     M2 --> C["<b>corrector</b><br/>Faraday, Ampère, Ohm<br/>→ B, E at n+1"]
     click M1 href "ppc_order.html#push-1-moments-only" "Push 1"
     click M2 href "ppc_order.html#push-2-for-real" "Push 2"

``SolverPPC::advanceLevel``, :src:`src/amr/solvers/solver_ppc.hpp#L378-L404`.

How it works
------------

The field stages
^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 18 42 40

   * - Stage
     - Computes
     - Code
   * - predictor 1
     - :math:`B^p` from :math:`B^n, E^n` (Faraday); :math:`J` from
       :math:`B^p` (Ampère); :math:`E^p` from :math:`B^p` and the ion moments
       at :math:`n` (Ohm)
     - ``predictor1_``, :src:`src/amr/solvers/solver_ppc.hpp#L409-L445`
   * - average
     - :math:`B^{avg} = (B^n + B^p)/2`, :math:`E^{avg} = (E^n + E^p)/2`; ghost
       values of :math:`E^{avg}` filled at :math:`t + \Delta t/2`
     - ``average_``, L534-L562
   * - predictor 2
     - :math:`B^p` from :math:`B^n, E^{avg}`; :math:`J`; :math:`E^p` with the
       moments predicted by push 1
     - ``predictor2_``, L448-L484
   * - corrector
     - :math:`B^{n+1}` from :math:`B^n, E^{avg}`; :math:`J`; :math:`E^{n+1}`
       with the moments of push 2
     - ``corrector_``, L488-L531

Both pushes use the averaged fields :math:`E^{avg}, B^{avg}`
(``moveIons_``, :src:`src/amr/solvers/solver_ppc.hpp#L604-L608`).

Push 1: moments only
^^^^^^^^^^^^^^^^^^^^

``IonUpdater::updateAndDepositDomain_``
(:src:`src/core/numerics/ion_updater/ion_updater.hpp#L184-L235`), for each
population:

#. copy the domain particles into a scratch array (L196);
#. push the copy, and deposit the copies that end in the patch's
   non-level-ghost box (its domain plus the cells of same-level neighbours)
   (L200-L205);
#. push a copy of the level ghost particles, and deposit those that enter the
   owned domain (L211-L233).

The real particle arrays are not touched: they stay at time :math:`n`. Then the
moments are completed (border sums, level ghost Old/New, totals:
:doc:`/dev/moments`). No particle changes patch.

Push 2: for real
^^^^^^^^^^^^^^^^

``IonUpdater::updateAndDepositAll_``
(:src:`src/core/numerics/ion_updater/ion_updater.hpp#L242-L311`), for each
population:

#. push the domain particles in place; those leaving the owned cells go to the
   end of the array (L258-L263);
#. leavers in a same-level neighbour's cells are copied to the patch ghost
   array; all leavers are removed from the domain array (L265-L289);
#. push the level ghost particles; those that enter owned cells are copied into
   the domain array, those that leave the ghost layer are erased (L291-L304);
#. deposit domain and patch ghost particles (L306-L309).

Then ``SolverPPC::moveIons_`` (:src:`src/amr/solvers/solver_ppc.hpp#L590-L666`)
completes the moments, lets neighbours grab the leavers
(:doc:`/dev/grabbing_particles`), computes the totals and fills the shared ghost
nodes.

.. mermaid::

   flowchart TD
     U["push and deposit<br/>(IonUpdater)"] --> B["sum density and flux<br/>over patch borders"]
     B --> L["add level ghost Old/New<br/>(time-interpolated)"]
     L --> G["push 2 only: neighbours<br/>grab the leavers;<br/>drop non-owned cells"]
     G --> T["totals: charge density,<br/>bulk velocity, charge flux"]
     T --> X["shared ghost nodes:<br/>keep the larger value"]
     click B href "moments.html" "Moments"
     click L href "level_ghost_particles.html" "Level ghost particles"
     click G href "grabbing_particles.html" "Grabbing particles"

The tricks
----------

Push 1 works on copies
^^^^^^^^^^^^^^^^^^^^^^

**Problem.** Push 1 only exists to predict the moments at :math:`n+1`. If it
moved the real particles, push 2 would start from the wrong positions.

**Fix.** ``updateAndDepositDomain_`` copies the domain and level ghost arrays
into ``tmp_particles_`` before pushing (L196, L212). The memory of that scratch
array is released when the hierarchy is regridded (``IonUpdater::reset``,
called from ``SolverPPC::onRegrid``,
:src:`src/amr/solvers/solver_ppc.hpp#L115-L119`).

Every particle gets both half pushes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** The Boris pusher moves the position in two half steps, with the
velocity update in between. If particles that leave the patch after the first
half step were set aside, a domain particle and its copy on the neighbour patch
would not receive the same number of half steps, and their positions would no
longer agree.

**Fix.** ``Boris::move`` does not select particles after the first half step
for domain particles (the first selector is a no-op, ``boxing.noop``); all of
them get both half steps, and the selection happens at the end
(:src:`src/core/numerics/pusher/boris.hpp#L94-L140`; comment at L103-L105).
Level ghost particles are an exception: those that leave the ghost box after
the first half step are not pushed further, since they will be erased anyway
(selector ``inGhostBox``,
:src:`src/core/numerics/ion_updater/ion_updater.hpp#L294-L296`).

**History.** :issue:`571` (Boris push desynchronization).

Particles change patch only in push 2
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** The exchange of particles between patches is a parallel
communication; doing it in push 1, whose particles are thrown away, would be
wasted.

**Fix.** ``fillIonGhostParticles`` and the single-owner clean-up run only when
the mode is not ``domain_only``
(:src:`src/amr/solvers/solver_ppc.hpp#L633-L656`). In push 1, leavers towards a
neighbour are still deposited (they are in the non-level-ghost box), so the
predicted moments are complete after the border sum.

Level ghost particles entering the level in push 1
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Copies of level ghost particles that enter the domain in push 1 are deposited
but not kept anywhere; the nodes at the level border get the rest of their
density from Old and New (comment at
:src:`src/core/numerics/ion_updater/ion_updater.hpp#L223-L230`;
:doc:`/dev/level_ghost_particles`).

Go deeper
---------

* :doc:`/dev/moments` — deposition and sums over patch borders.
* :doc:`/dev/grabbing_particles` — particles changing patch.
* :doc:`/dev/level_ghost_particles` — particles entering a refined level.

Related
-------

* :doc:`/dev/one_time_step` — where ``SolverPPC::advanceLevel`` sits.
* :doc:`/numerics/time_integration` — the scheme, for users.
* :doc:`/dev/electromagnetic` — Faraday, Ampère and Ohm.
* :doc:`/numerics/particles` — the Boris pusher and particle shapes.
