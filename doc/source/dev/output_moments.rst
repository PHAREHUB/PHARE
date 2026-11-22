.. _dev-output-moments:

================================
Ion moments computed for outputs
================================

.. container:: breadcrumb

   You are here: :doc:`/dev/index` › :doc:`/dev/diagnostics_restarts` › Ion moments for outputs

**In one sentence:** the time loop keeps only the ion moments Ohm's law reads; the
population fluxes on the moment nodes, their total, the bulk velocity and the momentum
tensors are computed by ``HybridLazyMoments`` when a diagnostic writes them or a Python
getter asks for them.

Why: since Ohm's law reads the ion current deposited on the E edges
(:doc:`/dev/moments`), these quantities only feed outputs. Computing them at every step
would cost time for nothing.

When they are computed
----------------------

.. plot::
   :include-source: False
   :caption: Five coarse steps. The time loop updates n and F_E at every step; the output
      moments are computed only when an output asks, once per level and time, and never
      when one level is written alone in the middle of a step.

   from phare_draw import *
   fig, (a,) = strip(1, ["who computes what, over five steps"], height=3.0, aspect=False)
   rows = {"time loop": 3, "diagnostic write": 2, "Python getter": 1, "dump_level": 0}
   for name, y in rows.items():
       a.plot([0.6, 5.4], [y, y], color=C["light"], lw=1)
       text(a, 0.5, y, name, ha="right", bold=True)
   for k in range(1, 6):
       mark(a, "primal", k, 3, color=C["rho"], size=9)
       text(a, k, -0.7, f"step {k}", color=C["grey"])
   text(a, 3, 3.4, "n, mass density, F_E: every step", color=C["rho"])
   for k in (2, 4):
       mark(a, "primal", k, 2, color=C["ok"], size=12)
   text(a, 3, 2.4, "compute, then write", color=C["ok"])
   mark(a, "primal", 3, 1, color=C["ok"], size=12)
   mark(a, "dual", 3.3, 1, color=C["grey"], size=10)
   text(a, 3.15, 1.4, "1st call computes, 2nd reuses", color=C["ok"])
   text(a, 4.5, 0, "✕", color=C["bad"], bold=True, fs=14)
   text(a, 4.5, 0.4, "mid-step: skipped", color=C["bad"])
   a.plot([-1.3, 5.6], [-1.0, 3.7], alpha=0)  # room for the row names
   finish(fig, [a], align=False)

* ``HybridLazyMoments`` is a member of the hybrid model, shared by the diagnostics and the
  Python data wrangler (:src:`src/amr/physical_models/hybrid_model.hpp#L60-L62`).
* Diagnostics: ``DiagnosticsManager::dump`` runs a diagnostic's ``compute`` whenever it
  is due to be computed or written
  (:src:`src/diagnostic/diagnostic_manager.hpp#L257-L280`). The H5 fluid writer asks for
  the flux if ``flux``, ``charge_flux`` or ``bulkVelocity`` is written, for the tensor if a
  ``momentum_tensor`` is (``FluidDiagnosticWriter::compute``,
  :src:`src/diagnostic/detail/types/fluid.hpp#L76-L100`); the VTK writer asks for the flux
  only (:src:`src/diagnostic/detail/vtk_types/fluid.hpp#L375-L393`). Both go through
  ``ModelView::fillLazyMoments``
  (:src:`src/diagnostic/diagnostic_model_view.hpp#L191-L197`).
* Python: the getters of population fluxes and bulk velocity call ``fillLazyFluxes_``
  first (:src:`src/python3/patch_level.hpp#L327-L335`). ``getChargeFluxE`` reads the E-edge
  ion current of the time loop directly
  (:src:`src/python3/patch_level.hpp#L260-L277`).
* Results are kept per level and per ions time: a second request for the same level object
  at the same time, for moments already computed, does nothing (``HybridLazyMoments::fill``,
  :src:`src/amr/physical_models/hybrid_lazy_moments.hpp#L66-L86`).

How one level is computed
-------------------------

.. flow::

   HybridLazyMoments::fill | nothing to do if this level was done at this time
   then
   each patch, each population | zero F (and M); deposit the domain particles and the level ghost particles Old the patch owns
   then
   border sums | add the neighbours' parts, in place, in one order; schedules built for this fill, then dropped
   then
   totals | charge flux on the moment nodes, bulk velocity, total momentum tensor

* The deposit uses the same kernels as the time loop; the population density it also
  produces goes to a scratch field, since the time loop owns the populations' densities
  (:src:`src/amr/physical_models/hybrid_lazy_moments.hpp#L94-L124`).
* Level ghost particles: each patch deposits only those it owns
  (``LevelGhostOwnership``, :doc:`/dev/level_ghost_particles`), from ``levelGhostParticlesOld``,
  which holds the particles at that time when every level is at the same time.
* Border sums: ``BorderSumRefiner`` with ``fillTransient``
  (:src:`src/amr/physical_models/hybrid_lazy_moments.hpp#L126-L132`); same in-place sums as
  the time loop (:doc:`/dev/moments`).
* Totals: ``computeChargeFluxOnMomentNodes``, ``computeBulkVelocity`` (divides by the mass
  density the time loop deposited), ``computeFullMomentumTensor``
  (:src:`src/amr/physical_models/hybrid_lazy_moments.hpp#L134-L143`).

The tricks
----------

Only when all levels are at the same time
   **Problem.** The level ghost particles deposited are ``levelGhostParticlesOld``, which hold
   the particles at the right time only when every level is at the same time.
   **Fix.** ``fill`` reads the ions' time stamp on every level and every rank, and throws on
   every rank if they differ; the per-level dumps (fine and emergency dumps) skip these
   quantities (:src:`src/amr/physical_models/hybrid_lazy_moments.hpp#L151-L175`,
   :doc:`/dev/diagnostics_restarts`).

Every write computes first
   **Problem.** A write due to elapsed wall-clock time, with no compute due, would write
   values the time loop no longer refreshes.
   **Fix.** ``dump`` computes a diagnostic whenever it writes it
   (:src:`src/diagnostic/diagnostic_manager.hpp#L263-L268`).

The getters are collective
   **Problem.** The border sums communicate between MPI ranks: a getter called on one rank
   only would wait forever.
   **Fix.** Population flux and bulk velocity getters must be called on every rank, as the
   docstrings say (:src:`pyphare/pyphare/data/wrangler.py#L20-L28`). The earliest and latest
   times are reduced over all ranks before any decision, a rank with no patch included, so
   they all take the same one (``mpi::min``,
   :src:`src/amr/physical_models/hybrid_lazy_moments.hpp#L166-L167`).

Schedules are not kept
   **Problem.** A kept schedule holds its level: after a regrid it would keep the replaced
   level alive, with its particles.
   **Fix.** The border-sum schedules are built for one fill and dropped; the ownership,
   made of boxes only, is rebuilt when the level object changes
   (:src:`src/amr/physical_models/hybrid_lazy_moments.hpp#L178-L220`).

Nothing new in restart files
   **Problem.** A new registered field is allocated on every patch and written to restart
   files.
   **Fix.** The population density goes into the messenger's scratch field
   ``PHARE_sumField``, which the VTK writer uses too; it is valid only inside ``fill``,
   which ``dump`` runs before any write
   (:src:`src/amr/physical_models/hybrid_lazy_moments.hpp#L241-L244`,
   :src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L1109-L1114`).

Charge flux on covered cells comes from the coarse particles
   The charge flux on the moment nodes is not coarsened at synchronization: on coarse cells
   under a finer level it comes from the coarse level's own particles, like the population
   fluxes and the bulk velocity. The charge density there is still coarsened
   (:src:`src/amr/messengers/hybrid_hybrid_messenger_strategy.hpp#L971-L986`).

Until computed, they hold old values
   Between outputs, these fields keep the values of the last computation, or NaN on a level
   that never had one (fields are allocated NaN)
   (:src:`src/core/data/ions/ions.hpp#L166-L171`). Read them only through a diagnostic or a
   getter.

Go deeper
---------

* :doc:`/dev/moments` — the deposit and the border sums this computation reuses.
* :doc:`/dev/level_ghost_particles` — which patch deposits each level ghost particle.

Related
-------

* :doc:`/dev/diagnostics_restarts` (up) · :doc:`/dev/python_to_cpp` ·
  :doc:`/user_guide/diagnostics`
