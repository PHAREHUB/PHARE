# Inner-Boundary Extraction Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Port the MHD inner-boundary (embedded body) feature from `embedded-boundary` onto the `inner-boundary` branch (based on master 3fcc2509), adapted to master's total-B formulation and periodic-only outer boundaries.

**Architecture:** Final-state copy + adapt. IB-owned files are copied at their `embedded-boundary` tip state (with the distance-aware-fill commit 967fa6d6 reverted), renamed to master's field names (B1→B, Etot1→Etot, B0 removed), and the integration hunks are re-fitted by hand into master's files. See `docs/superpowers/specs/2026-07-06-inner-boundary-extraction-design.md`.

**Tech Stack:** C++20 header templates, SAMRAI PatchData, pybind11/pyphare, CMake + CTest, gtest.

## Global Constraints

- Branch: `inner-boundary`; source of ported code: `embedded-boundary` (tip 9ca9e5bd).
- Build with `tools/cmake.sh` (NOT bare `uv run cmake --build`), `-j 12` max; existing `build/` dir is RelWithDebInfo.
- Run anything importing PHARE modules with `phare-run` alias semantics: venv python + `PYTHONPATH=$PHARE_HOME/build:$PHARE_HOME/pyphare:$PHARE_HOME`.
- Rename map applied to ALL ported code: `Etot1` → `Etot`; `B1` → `B` (identifier and `MHDQuantity::Vector::B1` → `Vector::B`, `statenew.B1` → `statenew.B`); every `B0`-specific slot/branch is deleted, not renamed.
- Nothing from the outer-BC lineage is ported: no `boundaryManager`, no `dict["grid"]["boundary_conditions"]`, no LODI/characteristic files.
- Distance-aware fill (967fa6d6) excluded: ghosts with non-interpolable mirror are skipped (pre-967fa6d6 behavior).
- After single-target builds, verify the target actually relinked (check mtime) before running tests — a failed ninja build leaves a stale binary.
- Commit messages end with `Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>`.

---

### Task 1: Dependency headers + core inner_boundary sources (part 1 of commit 1)

**Files:**
- Create (copy from `embedded-boundary:` same path):
  - `src/core/data/field/field_traits.hpp`
  - `src/core/data/grid/gridlayout_traits.hpp`
  - `src/core/data/tensorfield/tensorfield_traits.hpp`
  - `src/core/numerics/thermo/thermo.hpp`
  - `src/core/numerics/thermo/thermo_factory.hpp`
  - `src/core/numerics/primite_conservative_converter/conversion_utils.hpp`
  - `src/core/numerics/interpolator/field_at_point.hpp` (pre-967fa6d6 state)
  - `src/core/inner_boundary/*.hpp` — all 24 files
- Modify: `src/core/CMakeLists.txt` (add new headers to `SOURCES_INC`)

**Interfaces:**
- Produces: `core::InnerBoundaryManager<PhysicalQuantity, Field, GridLayout, VecField>::create(dict, scalarQuantities, vectorQuantities, thermo)` returning `std::unique_ptr` (null when no `inner_boundary` in dict); `core::Thermo` / `core::makeThermo`; `InnerBoundaryGeometry::characteristicLength()`; classifier + `InnerBoundaryMeshData` with degraded face/edge lists; `FieldAtPoint` interpolation.
- Consumes: master's `MHDQuantity` (`Scalar::rho`, `Scalar::Etot`, `Vector::B`, `Vector::E`, `Vector::rhoV`), `initializer/data_provider.hpp`, `layout.nbrGhosts()`.

- [ ] **Step 1: Verify which files 967fa6d6 must be reverted on, and that nothing later touched IB core**

```bash
git show --stat 967fa6d6            # lists the 9 IB files + field_at_point.hpp
git log --oneline 967fa6d6..embedded-boundary -- src/core/inner_boundary/ src/core/numerics/interpolator/field_at_point.hpp
```
Expected: second command prints NOTHING (no later commits). If it prints commits, inspect them and copy tip state instead, reverting 967fa6d6 with `git diff 967fa6d6 967fa6d6^ -- <file> | git apply` per file.

- [ ] **Step 2: Copy files.** For each file F touched by 967fa6d6, copy `git show 967fa6d6^:F`; for every other file, copy `git show embedded-boundary:F`:

```bash
mkdir -p src/core/inner_boundary src/core/numerics/thermo
for f in $(git ls-tree --name-only embedded-boundary src/core/inner_boundary/); do
  if git show 967fa6d6 --stat --name-only | grep -q "^$f$"; then src=967fa6d6^; else src=embedded-boundary; fi
  git show $src:$f > $f
done
git show 967fa6d6^:src/core/numerics/interpolator/field_at_point.hpp > src/core/numerics/interpolator/field_at_point.hpp
for f in src/core/data/field/field_traits.hpp src/core/data/grid/gridlayout_traits.hpp \
         src/core/data/tensorfield/tensorfield_traits.hpp src/core/numerics/thermo/thermo.hpp \
         src/core/numerics/thermo/thermo_factory.hpp \
         src/core/numerics/primite_conservative_converter/conversion_utils.hpp; do
  git show embedded-boundary:$f > $f
done
```

- [ ] **Step 3: Apply rename map** across `src/core/inner_boundary/` (see Global Constraints). Mechanical part:

```bash
grep -rn 'Etot1\|B1\|B0' src/core/inner_boundary/   # inventory first
sed -i 's/Etot1/Etot/g' src/core/inner_boundary/*.hpp
```
Then hand-edit the B1/B0 sites (4 files: `inner_boundary_manager.hpp`, `field_total_energy_from_pressure_inner_boundary_condition.hpp`, `field_ionospheric_convection_momentum_inner_boundary_condition.hpp`, `inner_boundary_condition_factory.hpp`):
- `inner_boundary_manager.hpp`: safe-state struct keeps ONE `std::array<double,3> B`; dict keys `B/x, B/y, B/z`; `safeVectors_[Vector::B] = s.B;` (delete the B0 line); total-energy safe value computed with `s.B` components.
- TEFP + ionospheric conditions: `ctx.statenew.B1` → `ctx.statenew.B`; face-to-cell-center projections unchanged (master B is face-centered).
- Delete any remaining `B0` logic (grep must return zero hits for `B0` afterwards).
- Doc comments mentioning "B1"/"Etot1"/"perturbation field" updated to total-B wording.
- In the ghost-fill loop where a non-interpolable mirror is skipped, keep/restore the `continue;` and add a comment: interpolation-dependent fills near patch edges are future work; magnetospheric cases currently circumvent this with a 0th-order Dirichlet condition.

- [ ] **Step 4: Fix conversion_utils/thermo integration.** Check `conversion_utils.hpp` and `thermo*.hpp` compile against master (they must not include B0/B1-era headers; apply the same rename map if they do). If `makeThermo` requires dict keys master doesn't populate (e.g. `eos`), make it default to ideal-gas with `heat_capacity_ratio` from `dict["to_primitive"]["heat_capacity_ratio"]` — the key master already populates.

- [ ] **Step 5: Register headers in `src/core/CMakeLists.txt`.** Copy the `inner_boundary/` block from the branch diff (`git diff master embedded-boundary -- src/core/CMakeLists.txt`), taking ONLY inner_boundary + traits + thermo + conversion_utils lines (skip `numerics/boundary_condition/field_*` outer-BC lines).

- [ ] **Step 6: Compile check.** Headers are template-only; real check happens with the unit tests in Task 2. Just confirm no stray includes:

```bash
grep -rn '#include' src/core/inner_boundary/ src/core/numerics/thermo/ | grep -v '<' | awk -F'"' '{print $2}' | sort -u | while read h; do [ -f "src/$h" ] || echo "MISSING src/$h"; done
```
Expected: no MISSING lines.

*(No commit yet — commit 1 lands at end of Task 2.)*

### Task 2: Unit tests port + build + commit 1

**Files:**
- Create: `tests/core/utilities/inner_boundary/CMakeLists.txt` + `test_*.cpp` (10 files from branch)
- Modify: whichever CMake file registers `tests/core/utilities/*` subdirectories (find with `grep -rn 'utilities/ghost_width_calculator\|utilities/box' tests/ res/cmake/ CMakeLists.txt`)

**Interfaces:**
- Consumes: Task 1 headers, master gtest setup, `add_no_mpi_phare_test` macro from `res/cmake/def.cmake`.

- [ ] **Step 1: Copy test sources and CMakeLists from `embedded-boundary:tests/core/utilities/inner_boundary/`** (all files; there are no separate distance-aware-fill test files, but individual EXPECTs may encode distance-aware behavior).

- [ ] **Step 2: Apply the rename map to the tests** (same sed + hand-edit for B1/B0/Etot1; safe-state dict keys become `B/x` etc.).

- [ ] **Step 3: Register the subdirectory** next to the existing `utilities` test registrations found in Step 0 grep.

- [ ] **Step 4: Configure + build the test targets**

```bash
tools/cmake.sh    # or: uv run cmake --build build -j 12 --target test-inner_boundary (then remaining test-field_* targets)
```
Expected: clean compile. Iterate on compile errors — most will be renamed-symbol or master-API drift (evalOnBox/ufuncing changes); fix in the ported headers, not by weakening tests.

- [ ] **Step 5: Run the suite**

```bash
cd build && uv run ctest -R 'inner_boundary' --output-on-failure
```
Expected: all PASS. Failures encoding distance-aware expectations (mirror off-patch cases) are changed to expect the skip behavior instead.

- [ ] **Step 6: Commit 1**

```bash
git add src/core/inner_boundary src/core/numerics/thermo src/core/numerics/interpolator/field_at_point.hpp \
        src/core/numerics/primite_conservative_converter/conversion_utils.hpp \
        src/core/data/field/field_traits.hpp src/core/data/grid/gridlayout_traits.hpp \
        src/core/data/tensorfield/tensorfield_traits.hpp src/core/CMakeLists.txt tests/core/utilities/inner_boundary <registration file>
git commit -m "feat(mhd): core inner-boundary infrastructure (geometry, classifier, field BCs, manager)"
```

### Task 3: AMR integration — PatchData, model, solver, degraded flux (commit 2)

**Files:**
- Create (copy + adapt from branch, same paths):
  - `src/amr/data/inner_boundary/ghost_elem_data.hpp`, `ghost_elem_data_factory.hpp`, `ghost_elem_variable.hpp`
  - `src/amr/resources_manager/ghost_elem_resource.hpp`
- Modify (extract IB-only hunks from `git diff master embedded-boundary -- <file>`):
  - `src/amr/resources_manager/resources_manager_utilities.hpp`
  - `src/amr/physical_models/mhd_model.hpp`
  - `src/amr/level_initializer/mhd_level_initializer.hpp`
  - `src/amr/solvers/solver_mhd.hpp`
  - `src/amr/solvers/solver_mhd_model_view.hpp`
  - `src/amr/solvers/time_integrator/compute_fluxes.hpp`
  - `src/amr/solvers/time_integrator/euler_using_computed_flux.hpp`
  - `src/core/numerics/godunov_fluxes/godunov_fluxes.hpp`
  - `src/core/numerics/finite_volume_euler/finite_volume_euler.hpp`
  - `src/core/numerics/constrained_transport/upwind_constrained_transport.hpp`
  - `src/core/numerics/primite_conservative_converter/to_primitive_converter.hpp`
  - `src/core/models/mhd_state.hpp` (likely NO change needed — verify its diff is all B0/B1)

**Interfaces:**
- Consumes: `InnerBoundaryManager` API from Task 1; master's templated `SolverMHD<MHDModel, AMR_Types, TimeIntegratorStrategy, Messenger>`.
- Produces: `MHDModel::hasInnerBoundary()`, `MHDModel::innerBoundaryManager`, `MHDModel::setupInnerBoundaryState(level, time)`; `Godunov::degrade_fluxes_near_inner_boundary(...)`; `ConstrainedTransport::degrade_E_near_inner_boundary(state, meshData)`; FVE fluid/cut-cell-only evolution when IB present.

**Hunk selection rule** (applies to every "Modify" file): from `git diff master embedded-boundary -- <file>`, take ONLY hunks containing these markers: `innerBoundary`, `ibm`, `hasInnerBoundary`, `GhostElem`, `ghost_elem`, `degrade_`, `setupInnerBoundaryState`, `inner_boundary`. SKIP hunks about: `B0`/`B1` splitting, `boundaryManager`/outer BC, reflux, adaptive `dt`/`computeStableDt`, tagging, `Etot1` renames outside IB hunks, diagnostics (`BTotal_`, `EtotTotal_`, `divB_diag_`, `tmpField_`, `tmpVec_`). Inside kept hunks apply the rename map; drop `model.state.B0x_Ez`/`B0y_Ez` arguments entirely (CT calls take `(state, ...)` per master's existing signature).

- [ ] **Step 1: Copy the 4 new AMR files, apply rename map.**

- [ ] **Step 2: `mhd_model.hpp`** — add manager member + creation in constructor:

```cpp
std::vector<core::MHDQuantity::Scalar> scalarQuantities
    = {core::MHDQuantity::Scalar::rho, core::MHDQuantity::Scalar::Etot};
std::vector<core::MHDQuantity::Vector> vectorQuantities = {
    core::MHDQuantity::Vector::B, core::MHDQuantity::Vector::E,
    core::MHDQuantity::Vector::rhoV};
innerBoundaryManager
    = inner_boundary_manager_type::create(dict, scalarQuantities, vectorQuantities, thermo);
if (innerBoundaryManager)
    resourcesManager->registerResources(*innerBoundaryManager);
```
plus `thermo` member (from Task 1 Step 4 decision), `allocate()` hunk, `hasInnerBoundary()`, `setupInnerBoundaryState()` (rename-mapped; no B0 masking). NO `boundaryManager`.

- [ ] **Step 3: remaining files per the hunk selection rule.** Work file by file; after each file `git diff` it and re-read for leftover B0/B1/outer-BC identifiers.

- [ ] **Step 4: Full build**

```bash
tools/cmake.sh   # full build, -j 12
```
Expected: green. This compiles every pybind permutation — long. Fix compile errors within the ported hunks.

- [ ] **Step 5: Run existing MHD + new IB test suites (regression + feature)**

```bash
cd build && uv run ctest -R 'mhd|inner_boundary' --output-on-failure -j 12
```
Expected: all PASS (master's MHD suites must be untouched by the port when no `inner_boundary` is configured).

- [ ] **Step 6: Commit 2**

```bash
git add src/amr src/core/numerics src/core/models
git commit -m "feat(mhd): integrate inner boundary into AMR model, solver, and numerics"
```

### Task 4: Python DSL + initializer plumbing (commit 3)

**Files:**
- Modify: `pyphare/pyphare/pharein/simulation.py`, `pyphare/pyphare/pharein/initialize/general.py`
- Test: `pyphare/pyphare_tests/pharein/simulation_test.py` (port IB test additions from branch)

**Interfaces:**
- Consumes: dict keys read by `InnerBoundaryManager::create` (Task 1): `simulation/inner_boundary/{name,shape,condition_type,density,pressure,center,radius,point,normal,inactive_safe_state/{density,pressure,velocity,B}}`.
- Produces: `ph.Simulation(inner_boundary={...})` kwarg.

- [ ] **Step 1: Extract the `inner_boundary` hunks** from `git diff master embedded-boundary -- pyphare/pyphare/pharein/simulation.py pyphare/pyphare/pharein/initialize/general.py` (skip outer-BC/`boundary_types`/tagging-halo/adaptive-dt hunks — `inner_boundary_no_refinement_halo` is NOT ported). In `general.py`, the safe-state vector loop becomes `for vec in ("velocity", "B"):`.

- [ ] **Step 2: Port the IB cases of `simulation_test.py`** from the branch (rename-mapped).

- [ ] **Step 3: Run**

```bash
cd build && uv run ctest -R 'py3.*simulation' --output-on-failure   # or the registered name via ctest -N | grep -i simulation
```
Expected: PASS.

- [ ] **Step 4: Commit 3**

```bash
git add pyphare
git commit -m "feat(pharein): inner_boundary configuration DSL"
```

### Task 5: Periodic functional smoke test (commit 4)

**Files:**
- Create: `tests/functional/mhd_inner_boundary/sphere_periodic.py`, `tests/functional/mhd_inner_boundary/CMakeLists.txt`
- Modify: the CMake file registering `tests/functional/*` subdirectories (find with `grep -rn 'functional' tests/CMakeLists.txt res/cmake/`)

**Interfaces:**
- Consumes: full stack from Tasks 1–4; `pyphare.pharesee.run.Run` for post-checks.

- [ ] **Step 1: Write the case**, starting from `git show embedded-boundary:tests/functional/mhd_inner_boundary/sphere_init.py`, adapted: periodic domain (default), MHD-only model options matching a master permutation from `res/sim/all.txt`, uniform initial flow (e.g. rho=1, V=(1,0,0), P=1, B=(0,0,1e-2)) around a centred sphere with a `dirichlet`-type (or the case's default) condition_type + `inactive_safe_state`, ~50 time steps, dump at final time. After `Simulator(...).run()`: assert via pharesee that rho/P on the finest level are finite and positive outside the body and that in-body cells hold the safe-state density. Wrap in `unittest.TestCase` like sibling functional MHD tests; reset `ph.global_vars.sim = None` in setUp.

- [ ] **Step 2: Register in CTest** with `phare_python3_exec` (or `phare_mpi_python3_exec` mirroring sibling MHD functional tests), guarded by `if(HighFive)` and exec level ≥ 11 like other heavy MHD cases; mirror `git show embedded-boundary:tests/functional/mhd_inner_boundary/CMakeLists.txt`.

- [ ] **Step 3: Run it directly first** (faster loop; kill leaked mpirun before each rerun: `pkill -f mpirun || true`):

```bash
env -C build PYTHONPATH=$PWD/build:$PWD/pyphare:$PWD .venv/bin/python ../tests/functional/mhd_inner_boundary/sphere_periodic.py
```
Expected: completes, assertions pass, no NaN/exception.

- [ ] **Step 4: Run via ctest** (`uv run ctest -R mhd_inner_boundary --output-on-failure`) if within configured `PHARE_EXEC_LEVEL_MAX`; otherwise note the guard in the commit message.

- [ ] **Step 5: Commit 4**

```bash
git add tests/functional/mhd_inner_boundary <registration file>
git commit -m "test(mhd): inner-boundary functional smoke case (sphere in periodic box)"
```

### Task 6: Final verification

- [ ] **Step 1: Full suite** `cd build && uv run ctest -j 12 --output-on-failure` — no regressions vs a master baseline (if a test fails, check it also fails on master before touching it).
- [ ] **Step 2: `git diff master...inner-boundary --stat`** — confirm no outer-BC, B0/B1, reflux, adaptive-dt, or tagging files slipped in.
- [ ] **Step 3: grep tree for leftovers** `grep -rn 'B0\|B1\|Etot1\|boundaryManager' src/core/inner_boundary src/amr/data/inner_boundary` → zero hits.
