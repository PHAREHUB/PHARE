.. _dev-ppc-order:

=====================================
Push and deposition order in the PPC
=====================================

.. container:: breadcrumb

   You are here: :doc:`/dev/index` › :doc:`/dev/one_time_step` › Push and
   deposition order

**In one sentence:** the hybrid solver ``SolverPPC`` advances one :term:`level`
by one step with a :term:`predictor-corrector scheme <PPC>` in which particles
are pushed twice: the first push only *predicts* the ion moments at the new
time, the second really moves the particles.

Why two pushes: :term:`Ohm's law` needs the moments at n+1, which need the
particles pushed with fields at n+½, which need E at n+1. The scheme predicts,
then corrects (physics: :doc:`/numerics/time_integration`).

One step, four stages
---------------------

.. plot::
   :include-source: False
   :caption: What is known after each stage. Filled: final value; open square:
      predicted value; dashed arrow: push of a throw-away copy.

   from phare_draw import *
   rows = {"B": 3, "E": 2, "moments": 1, "particles": 0}
   X = {"n": 0, "n+½": 2, "n+1": 4}
   def base(a):
       for name, y in rows.items():
           a.plot([0, 4], [y, y], color=C["light"], lw=1)
           text(a, -0.3, y, name, ha="right", bold=True)
       for lab, x in X.items():
           text(a, x, -0.6, lab, color=C["grey"])
       for y in rows.values():
           mark(a, "primal", 0, y, color=C["ink"], size=8)
       a.plot([-1.6, 5.6], [0, 4.2], alpha=0)
   fig, (a1, a2, a3, a4) = strip(4, [
       "Predictor 1: Faraday, Ampère, Ohm with moments at n",
       "Push 1, on a copy: predicted moments at n+1",
       "Predictor 2, then push 2: particles really move",
       "Corrector: B and E at n+1 with the moments of push 2"], height=1.9)
   base(a1)
   mark(a1, "dual", 4, 3, color=C["B"]); mark(a1, "dual", 4, 2, color=C["E"])
   mark(a1, "dual", 2, 3, color=C["B"]); mark(a1, "dual", 2, 2, color=C["E"])
   text(a1, 2, 3.4, "average", color=C["grey"])
   mark(a1, "primal", 0.2, 4.2, color=C["ink"], size=8); text(a1, 0.45, 4.2, "final", ha="left")
   mark(a1, "dual", 1.5, 4.2, color=C["ink"], size=8); text(a1, 1.75, 4.2, "predicted / temporary", ha="left")
   base(a2)
   mark(a2, "dual", 2, 3, color=C["B"]); mark(a2, "dual", 2, 2, color=C["E"])
   mark(a2, "dual", 4, 3, color=C["B"]); mark(a2, "dual", 4, 2, color=C["E"])
   arrow(a2, (0.2, 0), (3.8, 0), color=C["ion"], lw=1.5, ls="--", ms=10, text="copy pushed with E, B at n+½")
   mark(a2, "dual", 4, 1, color=C["rho"])
   text(a2, 4.4, 0, "copy dropped", ha="left", color=C["grey"])
   base(a3)
   mark(a3, "dual", 4, 3, color=C["B"]); mark(a3, "dual", 4, 2, color=C["E"])
   mark(a3, "dual", 2, 3, color=C["B"]); mark(a3, "dual", 2, 2, color=C["E"])
   arrow(a3, (0.2, 0), (3.8, 0), color=C["ion"], lw=2, ms=10, text="pushed, change patch")
   mark(a3, "primal", 4, 0, color=C["ion"]); mark(a3, "primal", 4, 1, color=C["rho"])
   base(a4)
   mark(a4, "primal", 4, 3, color=C["B"]); mark(a4, "primal", 4, 2, color=C["E"])
   mark(a4, "primal", 4, 1, color=C["rho"]); mark(a4, "primal", 4, 0, color=C["ion"])
   finish(fig, [a1, a2, a3, a4], gap=0.5)

* Predictor 1: B^p from B^n, E^n; J from B^p; E^p from Ohm with moments at n
  (``predictor1_``, :src:`src/amr/solvers/solver_ppc.hpp#L409-L445`).
* Average: E, B at n+½ = (n + predicted)/2; E ghosts filled at t+Δt/2
  (``average_``, :src:`src/amr/solvers/solver_ppc.hpp#L534-L562`).
* Push 1 / push 2: ``moveIons_`` in mode ``domain_only`` / ``all``
  (:src:`src/amr/solvers/solver_ppc.hpp#L590-L666`).
* Predictor 2 uses the predicted moments; corrector uses those of push 2
  (``predictor2_``, ``corrector_``, :src:`src/amr/solvers/solver_ppc.hpp#L448-L531`).
* After each push: border sums, level ghost Old/New, totals (:doc:`/dev/moments`);
  after push 2 only: neighbours grab the leavers (:doc:`/dev/grabbing_particles`).

The tricks
----------

Push 1 works on copies
^^^^^^^^^^^^^^^^^^^^^^

**Problem.** Moving the real particles in push 1 would start push 2 from the
wrong place.
**Fix.** ``updateAndDepositDomain_`` pushes a copy in ``tmp_particles_``
(:src:`src/core/numerics/ion_updater/ion_updater.hpp#L184-L235`); its memory is
released at regrid (``onRegrid``, :src:`src/amr/solvers/solver_ppc.hpp#L115-L119`).

Every particle gets both half pushes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** Setting leavers aside after the first half step would desynchronize
a particle from its copy on the neighbour patch.
**Fix.** Domain particles are not selected between the two half steps
(``boxing.noop``; :src:`src/core/numerics/pusher/boris.hpp#L94-L140`). :issue:`571`.

Particles change patch only in push 2
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** Exchanging particles in push 1, whose particles are thrown away,
would waste communication.
**Fix.** The grab and the single-owner clean-up run only in mode ``all``
(:src:`src/amr/solvers/solver_ppc.hpp#L633-L656`); moments are exchanged after
both pushes.

Level ghost entrants of push 1
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Copies of level ghost particles entering the domain in push 1 are deposited,
not kept (comment at :src:`src/core/numerics/ion_updater/ion_updater.hpp#L223-L230`).

Go deeper
---------

* :doc:`/dev/moments` · :doc:`/dev/grabbing_particles` · :doc:`/dev/level_ghost_particles`

Related
-------

* :doc:`/dev/one_time_step` · :doc:`/numerics/time_integration` ·
  :doc:`/dev/electromagnetic` · :doc:`/numerics/particles`
