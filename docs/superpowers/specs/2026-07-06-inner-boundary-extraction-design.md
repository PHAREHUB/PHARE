# Inner-boundary feature extraction onto master

**Date:** 2026-07-06
**Branch:** `inner-boundary`, based on `master` (3fcc2509)
**Source:** `embedded-boundary` branch tip (9ca9e5bd)

## Goal

Isolate the MHD inner-boundary (embedded body) feature developed on `embedded-boundary`
onto a clean, PR-ready branch based on current master. Master lacks two features the
source branch was built on — magnetic-field splitting (B = B0 + B1) and outer physical
boundary conditions — so the ported code is adapted to master's total-B formulation and
periodic-only outer boundaries.

## Extraction mechanism

Final-state copy + adapt. The `embedded-boundary` history (161 commits over merge base,
merge-heavy) makes cherry-picking or history rewriting impractical. Instead:

1. Copy the IB-owned directories at their final state on `embedded-boundary`.
2. Revert the distance-aware-fill change (967fa6d6) on those files (see Omissions).
3. Hand-port the integration hunks into master's files, adapting names and APIs.

## Components ported

- `src/core/inner_boundary/` — geometry (sphere, plane), mesh classifier, mesh data,
  ghost-elem pack, `InnerBCContext`, field BC types (none, Dirichlet, Neumann,
  symmetric, antisymmetric, adaptive Dirichlet-or-Neumann, total-energy-from-pressure,
  ionospheric-convection momentum), condition factories, `InnerBoundaryManager` with
  config-driven `inactive_safe_state` (`setSafeState`).
- `src/amr/data/inner_boundary/` — per-patch ghost-elem SAMRAI PatchData (variable,
  factory, data) + resources-manager glue (`ghost_elem_resource.hpp`,
  `resources_manager_utilities.hpp` hunks).
- Degraded first-order ideal flux near under-resolved inner boundary:
  `InnerBoundaryGeometry::characteristicLength()`, classifier under-resolution flag +
  1-cell shell + degraded face/edge lists, `Godunov::degrade_fluxes_near_inner_boundary`,
  CT counterpart, two-pass solver hooks.
- Integration hunks (IB-only, re-fitted to master):
  `src/amr/physical_models/mhd_model.hpp` (manager ownership, allocate/register,
  `setupInnerBoundaryState`), `src/amr/solvers/solver_mhd.hpp`,
  `src/amr/solvers/solver_mhd_model_view.hpp` (BC application points),
  `src/amr/solvers/time_integrator/compute_fluxes.hpp`,
  `src/amr/solvers/time_integrator/euler_using_computed_flux.hpp`,
  `src/amr/level_initializer/mhd_level_initializer.hpp`,
  `src/core/numerics/godunov_fluxes/godunov_fluxes.hpp`,
  `src/core/numerics/finite_volume_euler/finite_volume_euler.hpp`,
  `src/core/numerics/constrained_transport/upwind_constrained_transport.hpp`,
  `src/core/numerics/primite_conservative_converter/to_primitive_converter.hpp`,
  `src/core/models/mhd_state.hpp`, `src/core/CMakeLists.txt`.
- Python DSL: `pyphare/pyphare/pharein/simulation.py` (`inner_boundary` kwarg
  validation), `pyphare/pyphare/pharein/initialize/general.py` (populateDict:
  name/shape/condition_type/density/pressure/center/radius/point/normal +
  `inactive_safe_state`).
- C++ initializer dict reads for the same keys.

## Adaptations

- **Total-B formulation:** `B1 -> B`, `Etot1 -> Etot` throughout; safe state carries a
  single `B` vector (DSL drops the `B0`/`B1` keys); all B0-specific logic removed
  (B0 masking, B0 safe-state slots, B0 face projections). Master's MHD state
  (`rho, V, B, P, rhoV, Etot, J, E`) is the direct target.
- **Distance-aware fill omitted:** the off-patch-mirror lever-arm extrapolation
  (967fa6d6) is not correct at this time. Restore the prior behavior: ghosts whose
  mirror point is not interpolable are skipped (the sibling-component
  interpolability filter is kept). A code comment records the current workaround for
  magnetospheric cases: 0th-order Dirichlet where fluid-value interpolation would be
  required. Fixing interpolation-dependent fills is future work.
- **Tagging refinement halo omitted:** `inner_boundary_no_refinement_halo` and the
  `concrete_tagger.hpp` changes stay behind; master's tagger is untouched. The
  ghost-classification fixes bundled in the same source commits are kept.
- **Solver shape:** hooks are written against master's templated `SolverMHD`
  (`TimeIntegratorStrategy` template parameter); the runtime-dispatch refactor from
  `remove-mhd-template-params` is not ported.
- **No outer-BC coupling:** nothing from the `boundary_conditions` lineage is ported;
  domains remain periodic.

## Tests

- Unit: `tests/core/utilities/inner_boundary/` full suite, minus distance-aware-fill
  tests, with B0/B1 assertions adapted to total B.
- Functional smoke: sphere inner boundary in a fully periodic box with an initial
  uniform flow, registered behind the heavy-test exec-level guard. Asserts: run
  completes without NaN/exception, inactive (in-body) cells hold the configured safe
  state, divB stays bounded.

## Commit plan

1. core IB infrastructure + unit tests
2. AMR integration (PatchData, model, solver, level initializer, degraded flux)
3. Python DSL + initializer plumbing
4. functional smoke test

## Risks

- Master's five post-merge-base commits (ufuncing #1192, patch-data-transfer #1201,
  evalOnBox #1200, mhd #1218) changed the very files the hunks land in — each hunk is
  re-fitted, not blindly copied.
- Ghost-width machinery differs from the source branch; the classifier sizes its shell
  from `layout.nbrGhosts()`, which exists on master, so this should carry over — verify
  against master's MHD ghost width at build time.
