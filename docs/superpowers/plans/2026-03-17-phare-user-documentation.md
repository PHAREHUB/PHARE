# PHARE End-User Documentation Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fill and expand PHARE's Sphinx documentation from mostly-empty stubs into comprehensive end-user docs covering theory, building, configuration, tutorials, and data analysis.

**Architecture:** All documentation lives in `doc/source/` as RST files. The existing Sphinx infrastructure (conf.py, Makefile, RTD theme) is unchanged. We fill empty pages, expand partial ones, create new tutorial and analysis pages, and update `index.rst` to add the tutorials toctree section.

**Tech Stack:** Sphinx RST, MathJax for equations, `.. code-block:: python` for code examples, `.. autoclass::` for API autodoc where appropriate.

**Spec:** `docs/superpowers/specs/2026-03-17-phare-user-documentation-design.md`

---

## File Structure

All paths relative to `doc/source/`:

**Existing files to expand/rewrite:**
| File | Current state | Action |
|------|--------------|--------|
| `index.rst` | Has toctrees for Theory/Building/Usage/Analysis | Add Tutorials toctree, remove `usage/examples` |
| `getting.rst` | Just git clone command | Expand with prerequisites, dependencies |
| `build.rst` | Two build commands | Expand with CMake options table, PYTHONPATH, troubleshooting |
| `usage/simulation_inputs.rst` | Partial (has Simulation, MaxwellianFluidModel, diagnostics basics) | Rewrite: complete parameter reference, add MHD, LoadBalancer |
| `usage/run_from_python.rst` | Basic (has PYTHONPATH, simple example) | Fix typos, add MPI execution, lifecycle, callbacks |
| `usage/run_from_exe.rst` | Empty (title only) | Stub redirecting to run_from_python |
| `theory/spatial_discretization.rst` | Empty (title only) | Fill with Yee lattice, stencils, ghost cells |
| `theory/temporal_discretization.rst` | Empty (title only) | Fill with predictor-corrector, Boris pusher, subcycling |
| `theory/amr.rst` | Empty section headers | Fill all sections with content |
| `pharesee/get_data.rst` | Partial (Run basics, incomplete hierarchy) | Rewrite with complete API reference |
| `pharesee/plotting_fields.rst` | Empty (title only) | Fill with 1D/2D plotting examples |
| `pharesee/plotting_distributions.rst` | Empty (title only) | Fill with phase space, velocity distributions |

**New files to create:**
| File | Content |
|------|---------|
| `theory/mhd.rst` | MHD equations, reconstruction, Riemann solvers, time stepping |
| `tutorials/alfven_wave_1d.rst` | Beginner: 1D Alfven wave |
| `tutorials/harris_2d.rst` | Intermediate: 2D reconnection |
| `tutorials/multi_population_1d.rst` | Intermediate: ion-ion beam |
| `tutorials/mhd_shock_1d.rst` | MHD beginner: 1D shock tube |
| `tutorials/mhd_orszag_tang_2d.rst` | MHD intermediate: 2D vortex |
| `tutorials/amr_tagging.rst` | Advanced: adaptive refinement |
| `pharesee/reconnection_analysis.rst` | Reconnection-specific analysis tools |

**Files to remove:**
| File | Reason |
|------|--------|
| `usage/examples.rst` | Replaced by `tutorials/` section |
| `usage/maxwellian_fluid_model.rst` | Content absorbed into `usage/simulation_inputs.rst` (was an orphan autodoc page) |

**Reference files** (read for accuracy, do not modify):
- `pyphare/pyphare/pharein/simulation.py` — all Simulation parameters and defaults
- `pyphare/pyphare/pharein/maxwellian_fluid_model.py` — MaxwellianFluidModel API
- `pyphare/pyphare/pharein/mhd_model.py` — MHDModel API
- `pyphare/pyphare/pharein/electron_model.py` — ElectronModel API
- `pyphare/pyphare/pharein/diagnostics.py` — all diagnostics classes and valid quantities
- `pyphare/pyphare/pharein/load_balancer.py` — LoadBalancer parameters
- `pyphare/pyphare/pharesee/run/run.py` — Run class methods and signatures
- `pyphare/pyphare/simulator/simulator.py` — Simulator class
- `src/core/numerics/ohm/ohm.hpp` — Ohm's law implementation
- `src/core/numerics/faraday/faraday.hpp` — Faraday implementation
- `src/core/numerics/ampere/ampere.hpp` — Ampere implementation
- `src/core/numerics/pusher/boris.hpp` — Boris pusher algorithm
- `src/core/data/grid/gridlayout_impl_yee.hpp` — Yee grid centering
- `src/amr/solvers/solver_ppc.hpp` — predictor-corrector algorithm
- `src/amr/multiphysics_integrator.hpp` — subcycling logic
- `tests/functional/alfven_wave/alfven_wave1d.py` — tutorial 1 source
- `tests/functional/harris/harris_2d.py` — tutorial 2 source
- `tests/functional/ionIonBeam/ion_ion_beam1d.py` — tutorial 3 source
- `tests/functional/mhd_shock/mhd_shock.py` — tutorial 4 source
- `tests/functional/mhd_orszagtang/orszag_tang.py` — tutorial 5 source
- `tests/functional/translation/translat1d.py` — tutorial 6 source

---

## Task 1: Update index.rst and create tutorials directory

**Files:**
- Modify: `doc/source/index.rst`
- Delete: `doc/source/usage/examples.rst`

- [ ] **Step 1: Read current index.rst**

Read `doc/source/index.rst` to understand the current toctree structure.

- [ ] **Step 2: Update index.rst toctree**

Replace the content of `doc/source/index.rst`. Changes:
1. Add `theory/mhd` to the Theory toctree
2. Remove `usage/examples` from the Usage toctree
3. Add a new TUTORIALS toctree section between USAGE and ANALYSIS with entries: `tutorials/alfven_wave_1d`, `tutorials/harris_2d`, `tutorials/multi_population_1d`, `tutorials/mhd_shock_1d`, `tutorials/mhd_orszag_tang_2d`, `tutorials/amr_tagging`
4. Add `pharesee/reconnection_analysis` to the Analysis toctree
5. Remove the "Warning: This documentation is a work in progress" line

- [ ] **Step 3: Create tutorials directory**

```bash
mkdir -p doc/source/tutorials
```

- [ ] **Step 4: Delete obsolete files**

```bash
rm doc/source/usage/examples.rst
rm doc/source/usage/maxwellian_fluid_model.rst
```

- [ ] **Step 5: Commit**

```bash
git add doc/source/index.rst doc/source/tutorials/
git rm doc/source/usage/examples.rst doc/source/usage/maxwellian_fluid_model.rst
git commit -m "docs: update index toctree, add tutorials section, remove obsolete pages"
```

---

## Task 2: Expand getting.rst and build.rst

**Files:**
- Modify: `doc/source/getting.rst`
- Modify: `doc/source/build.rst`

- [ ] **Step 1: Read current getting.rst and build.rst**

Read both files to see existing content.

- [ ] **Step 2: Rewrite getting.rst**

Replace content with:
- Title: "Getting PHARE"
- **Prerequisites** section with a table listing: C++20 compiler (gcc >= 10 or clang >= 13), CMake >= 3.20.1, MPI (OpenMPI/MPICH/Intel MPI), parallel HDF5 (built with `--enable-parallel`), Python >= 3.8. Use RST `.. list-table::` format.
- **Auto-downloaded dependencies** section: note that SAMRAI, pybind11, HighFive, GoogleTest are fetched automatically during build — users do NOT need to install them.
- **Getting the source** section: `git clone --recursive https://github.com/PHAREHUB/PHARE`
- **Python dependencies** section:
  ```bash
  python3 -m venv phare_venv
  source phare_venv/bin/activate
  pip install -r requirements.txt
  ```
  List key packages: numpy, scipy, h5py, matplotlib, dill, pyyaml

- [ ] **Step 3: Rewrite build.rst**

Replace content with:
- Title: "Building PHARE"
- **Production build** section with optimized cmake + make commands
- **Debug build** section with debug flags
- **CMake options** section with a table (Name, Default, Description) covering all flags from the spec: `devMode`, `test`, `testMPI`, `PHARE_MPI_PROCS`, `lowResourceTests`, `asan`, `ubsan`, `coverage`, `bench`, `SAMRAI_ROOT`, `PHARE_EXEC_LEVEL_MIN`, `PHARE_EXEC_LEVEL_MAX`
- **Setting up the environment** section with PYTHONPATH export:
  ```bash
  export PYTHONPATH=/path/to/PHARE/pyphare:/path/to/build:$PYTHONPATH
  ```
  Emphasize this is required before running any PHARE script.
- **Verifying the build** section with a minimal smoke test:
  ```python
  python3 -c "import pyphare.pharein as ph; print('PHARE ready')"
  ```
- **Troubleshooting** section with subsections for common issues:
  - HDF5 not parallel: check with `h5cc -showconfig | grep "Parallel HDF5"`
  - Wrong Python: use `-DPython_EXECUTABLE=/path/to/python3`
  - MPI not found: set `MPI_HOME`
  - SAMRAI download fails behind proxy: use `-DSAMRAI_ROOT=/path/to/prebuilt`

- [ ] **Step 4: Commit**

```bash
git add doc/source/getting.rst doc/source/build.rst
git commit -m "docs: expand getting started and build instructions with prerequisites, CMake options, troubleshooting"
```

---

## Task 3: Expand usage/run_from_python.rst and stub run_from_exe.rst

**Files:**
- Modify: `doc/source/usage/run_from_python.rst`
- Modify: `doc/source/usage/run_from_exe.rst`

- [ ] **Step 1: Read current run_from_python.rst and Simulator source**

Read `doc/source/usage/run_from_python.rst` and `pyphare/pyphare/simulator/simulator.py` for accurate API details.

- [ ] **Step 2: Rewrite run_from_python.rst**

Keep the existing structure but expand. Key additions:
- Fix typo "minimmum" → "minimum"
- **Running with MPI** section:
  ```bash
  mpirun -np 4 python3 -Ou my_script.py
  ```
  Explain `-Ou` enables unbuffered stdout/stderr for cleaner MPI output.
- **The Simulator lifecycle** section explaining the phases:
  1. `Simulator(sim)` — bind to Simulation config
  2. `sim.initialize()` — SAMRAI creates mesh hierarchy, loads initial conditions, dumps first diagnostics
  3. `sim.advance(dt)` — advance one timestep (or `sim.run()` to run until `final_time`)
  4. `sim.reset()` — clean up C++ objects
  Show a complete example script demonstrating this.
- **Post-advance callbacks** section:
  ```python
  def my_callback(time):
      print(f"Advanced to t={time}")
  simulator = Simulator(sim, post_advance=my_callback)
  ```
- **Dry run mode** section: `dry_run=True` in Simulation or `PHARE_DRY_RUN=1` environment variable. Validates config without running.
- **Log files** section: when `log_to_file=True` (default), each MPI rank writes to a separate log file.
- Keep the `.. autoclass::` directive for Simulator

- [ ] **Step 3: Rewrite run_from_exe.rst**

Replace with:
- Title: "Running from a compiled binary"
- A note that PHARE is designed to be run from Python:
  ```rst
  .. note::
     PHARE does not provide a standalone executable. All simulations are configured
     and launched from Python scripts. See :doc:`run_from_python` for details.
  ```

- [ ] **Step 4: Commit**

```bash
git add doc/source/usage/run_from_python.rst doc/source/usage/run_from_exe.rst
git commit -m "docs: expand run_from_python with MPI, lifecycle, callbacks; stub run_from_exe"
```

---

## Task 4: Rewrite usage/simulation_inputs.rst

This is the largest single page. It becomes the comprehensive parameter reference.

**Files:**
- Modify: `doc/source/usage/simulation_inputs.rst`

- [ ] **Step 1: Read source files for accurate defaults**

Read these files to get exact parameter names, types, and defaults:
- `pyphare/pyphare/pharein/simulation.py` — Simulation class
- `pyphare/pyphare/pharein/maxwellian_fluid_model.py` — MaxwellianFluidModel
- `pyphare/pyphare/pharein/mhd_model.py` — MHDModel
- `pyphare/pyphare/pharein/electron_model.py` — ElectronModel
- `pyphare/pyphare/pharein/diagnostics.py` — all diagnostics classes
- `pyphare/pyphare/pharein/load_balancer.py` — LoadBalancer

- [ ] **Step 2: Write the Simulation block section**

**Important accuracy notes for this step:**
- Use `"modified_boris"` (snake_case) for `particle_pusher` default — the class docstring says `"modifiedBoris"` but the actual validation code uses `"modified_boris"`
- `hyper_resistivity` default is `0.0001` (non-zero!) — highlight this
- `tagging_threshold` default is `0.1`
- `nesting_buffer` default is `0` (no gap between levels)
- `clustering` default is `"tile"`

Rewrite the opening and Simulation block. Structure as:

1. Brief intro explaining the Python configuration approach
2. **The Simulation block** with parameter tables organized by category. Each table has columns: Name | Type | Default | Description. Categories:
   - Domain & mesh (`cells`, `dl`, `domain_size`, `boundary_types`)
   - Time stepping (`time_step`, `time_step_nbr`, `final_time`) with a note about the 2-of-3 rule
   - Numerics (`interp_order`, `particle_pusher`, `layout`, `refined_particle_nbr`)
   - AMR (`max_nbr_levels`, `refinement`, `refinement_boxes`, `tagging_threshold`, `smallest_patch_size`, `largest_patch_size`, `nesting_buffer`, `tag_buffer`, `clustering`)
   - Physics (`resistivity`, `hyper_resistivity`, `hyper_mode`)
   - MHD-specific (`model_options`, `gamma`, `eta`, `nu`, `hall`, `res`, `hyper_res`, `reconstruction`, `limiter`, `riemann`, `mhd_timestepper`, `max_mhd_level`)
   - Output (`diag_options`, `restart_options`)
   - Misc (`strict`, `dry_run`, `description`, `write_reports`)

Include a complete example Simulation block after the tables.

Cross-reference: link to :doc:`../theory/amr` for AMR details.

- [ ] **Step 3: Write MaxwellianFluidModel section**

Keep the existing Harris sheet examples (they're good). Add:
- Explanation that `bx`, `by`, `bz` are callables of spatial coordinates
- Table of per-population parameters: `charge`, `mass`, `density`, `vbulkx/y/z`, `vthx/y/z`, `nbr_part_per_cell`, `init` (dict with `seed`), `density_cut_off`
- Note about function signatures per dimension: `f(x)` for 1D, `f(x, y)` for 2D, `f(x, y, z)` for 3D
- Warning that functions must be periodic (PHARE validates at boundaries)
- Multi-population example (keep the existing one)

- [ ] **Step 4: Write MHDModel section**

New subsection covering:
- Parameters as spatial callables. Defaults: `density=1.0`, `vx=1.0`, `vy=0.0`, `vz=0.0`, `bx=1.0`, `by=0.0`, `bz=0.0`, `p=1.0`. Note: `vx`, `density`, `bx`, `p` default to 1.0 while others default to 0.0 — users should always set these explicitly.
- Must set `model_options=["MHDModel"]` in Simulation
- Guidance box on choosing reconstruction/limiter/Riemann solver:
  - Start with WENOZ + Rusanov for robustness
  - HLLD for less diffusion once stable
  - MP5 for highest accuracy on smooth problems
- Complete MHD example with all parameters

- [ ] **Step 5: Write ElectronModel section**

Expand existing brief section:
- `closure="isothermal"` (only option currently)
- `Te` parameter: electron temperature in simulation units
- Physical meaning: Pe = n * Te everywhere, where n is the plasma density

- [ ] **Step 6: Write Diagnostics section**

Rewrite the diagnostics section with:
- Common parameters table: `write_timestamps`, `flush_every` (default 1 — flushes every dump, safe but costly for frequent writes), `quantity`
- Note that `write_timestamps` is the only reliably working timestamp mechanism. `elapsed_timestamps` exists in the code but has a validation bug that may reject it — do not document it as a supported interface.
- **ElectromagDiagnostics**: quantity "E" or "B"
- **FluidDiagnostics** with two sub-tables:
  - Total ion quantities (no `population_name`): "charge_density", "mass_density", "bulkVelocity", "momentum_tensor"
  - Per-population (requires `population_name`): "flux", "momentum_tensor", "density"
  - Note about "pressure_tensor" as a convenience wrapper
- **ParticleDiagnostics**: "domain", "levelGhost", "space_box" with `extent`. Warning about output size.
- **MHDDiagnostics**: "rho", "V", "P", "rhoV", "Etot"
- **MetaDiagnostics**: "tags"
- **InfoDiagnostics**: "particle_count"
- Diagnostic output formats: "phareh5" (default) and "pharevtkhdf"
- Timestamp generation patterns with code examples

- [ ] **Step 7: Write LoadBalancer section**

New subsection:
- Parameter table: `active`, `mode`, `tol`, `auto`, `every`, `on_init`
- Expert parameters: `next_rebalance`, `next_rebalance_backoff_multiplier`, `max_next_rebalance`
- Brief guidance: use nppc mode for PIC, homogeneous for MHD
- Example

- [ ] **Step 8: Commit**

```bash
git add doc/source/usage/simulation_inputs.rst
git commit -m "docs: comprehensive rewrite of simulation_inputs with complete parameter reference, MHD, LoadBalancer"
```

---

## Task 5: Fill theory/spatial_discretization.rst

**Files:**
- Modify: `doc/source/theory/spatial_discretization.rst`

- [ ] **Step 1: Read Yee grid implementation for accuracy**

Read `src/core/data/grid/gridlayout_impl_yee.hpp` to understand exact centering for each field component and each dimension.

- [ ] **Step 2: Write spatial_discretization.rst**

Replace the stub with full content:

- **The Yee lattice** — explain staggered grid: fields are not all at the same location within a cell. Each component of E, B, J has a specific centering (primal or dual) in each direction.
- **Field centering table** — RST table showing centering (primal/dual) for each scalar quantity (Bx, By, Bz, Ex, Ey, Ez, Jx, Jy, Jz, rho) in each direction (x, y, z). Derive from `GridLayoutImplYee`.
- **Finite difference operators** — show the curl stencils used for Faraday (∂B/∂t = -∇×E) and Ampere (μ₀J = ∇×B). Write out the 1D and 2D forms explicitly with equations using `.. math::` directives.
- **Gradient operator** — used in electron pressure term (∇Pe/n). Show the centered difference formula.
- **Ghost cells** — explain their purpose (stencil support at patch boundaries). Width depends on interpolation order: order 1 needs 1 ghost, order 2 needs 2, order 3 needs 3. Include a diagram (ASCII art in RST) showing physical cells vs ghost cells in 1D.
- Cross-references: link to :doc:`hybridpic` for the equations, :doc:`../usage/simulation_inputs` for `interp_order` parameter.

- [ ] **Step 3: Commit**

```bash
git add doc/source/theory/spatial_discretization.rst
git commit -m "docs: fill spatial discretization theory - Yee lattice, stencils, ghost cells"
```

---

## Task 6: Fill theory/temporal_discretization.rst

**Files:**
- Modify: `doc/source/theory/temporal_discretization.rst`

- [ ] **Step 1: Read solver and pusher implementations**

Read `src/amr/solvers/solver_ppc.hpp` (advanceLevel method) and `src/core/numerics/pusher/boris.hpp` for algorithm details. Also read `src/amr/multiphysics_integrator.hpp` for subcycling logic (look for `ratio.max() * ratio.max()`).

- [ ] **Step 2: Write temporal_discretization.rst**

Replace the stub with:

- **The predictor-corrector scheme** — overview paragraph explaining that PHARE uses a second-order predictor-corrector to advance fields and particles self-consistently.
- **Algorithm steps** — numbered list with equations:
  1. Predictor 1: Faraday B→B*, Ampere B*→J*, Ohm (n,Ve,Pe,B*,J*)→E*
  2. Average: B̄=(B+B*)/2, Ē=(E+E*)/2
  3. Move ions: push domain particles with B̄, Ē; deposit moments
  4. Predictor 2: repeat predictor with updated moments
  5. Average: update B̄, Ē
  6. Move ions: push all particles (domain + ghosts)
  7. Corrector: Faraday→B_final, Ampere→J_final, Ohm→E_final
- **The Boris pusher** — explain the 3-step Boris algorithm for advancing particle velocity under Lorentz force:
  1. Half-step electric acceleration: v⁻ = vⁿ + (q/m)(E)(dt/2)
  2. Magnetic rotation: t = (q/m)(B)(dt/2), s = 2t/(1+t²), v' = v⁻ + v⁻×t, v⁺ = v⁻ + v'×s
  3. Half-step electric acceleration: vⁿ⁺¹ = v⁺ + (q/m)(E)(dt/2)
  4. Position update: xⁿ⁺¹ = xⁿ + vⁿ⁺¹ dt
  Use `.. math::` for the equations.
- **Stability constraint** — particles must not cross more than one cell per timestep. If they do, the code throws an error. Users must choose `time_step` small enough.
- **AMR subcycling** — finer levels take more, smaller timesteps. The time step ratio is the *square* of the refinement ratio: dt_fine = dt_coarse / ratio². With ratio=2, each fine level takes 4 substeps per coarse step. This is because whistler waves impose dt ~ dx² scaling. Link to :doc:`amr`.

- [ ] **Step 3: Commit**

```bash
git add doc/source/theory/temporal_discretization.rst
git commit -m "docs: fill temporal discretization - predictor-corrector, Boris pusher, subcycling"
```

---

## Task 7: Fill theory/amr.rst

**Files:**
- Modify: `doc/source/theory/amr.rst`

- [ ] **Step 1: Read AMR-related source files**

Read `src/amr/multiphysics_integrator.hpp`, `src/amr/messengers/hybrid_messenger_strategy.hpp`, and `src/amr/resources_manager/resources_manager.hpp` for implementation details.

- [ ] **Step 2: Write amr.rst**

Keep the existing title and section headers but fill all sections:

- **Introduction** paragraph: what AMR is, why it matters for hybrid PIC (resolving sub-ion structures like current sheets without paying the cost everywhere).
- **Patch based approach**: SAMRAI library provides the AMR framework. The domain is divided into rectangular patches. Multiple refinement levels form a hierarchy: Level 0 (coarsest) covers the whole domain, Level 1+ cover subregions at higher resolution. Refinement ratio is fixed at 2 (cell size halves at each level).
- **Recursive time integration**: SAMRAI drives recursive subcycling. Coarsest level advances one step, then each finer level advances ratio² = 4 substeps. Synchronization happens bottom-up after subcycling. Link to :doc:`temporal_discretization`.
- **Two refinement strategies**:
  - Static boxes: user specifies `refinement_boxes` dict mapping level → list of Box regions
  - Adaptive tagging: `refinement="tagging"`, code evaluates B-field gradients, cells exceeding `tagging_threshold` are tagged for refinement
- **Field refinement**: when a finer level is created, field values are interpolated from the coarser level using spatial interpolation. For subcycling, time interpolation between coarse-level snapshots provides ghost data at intermediate times.
- **Particle refinement**: `refined_particle_nbr` controls how many fine particles replace one coarse particle at refinement boundaries. Particles are split with adjusted weights to conserve moments.
- **Field coarsening**: after fine-level advancement, fine data is coarsened (averaged) back to the coarse level to keep levels consistent.
- **Fields at level boundaries**: ghost cells at patch/level boundaries are filled from neighboring patches (same level) or interpolated from the coarser level.
- **Particles at level boundaries**: level ghost particles are created from coarser level via splitting. Patch ghost particles come from neighboring patches at the same level.
- **Patch size constraints**: `smallest_patch_size` (minimum cells per patch, default depends on `interp_order`), `largest_patch_size` (maximum), `nesting_buffer` (gap between levels). `clustering` algorithm: "tile" (default) or "berger" (Berger-Rigoutsos).

- [ ] **Step 3: Commit**

```bash
git add doc/source/theory/amr.rst
git commit -m "docs: fill AMR theory - patch hierarchy, subcycling, refinement strategies, particle splitting"
```

---

## Task 8: Create theory/mhd.rst

**Files:**
- Create: `doc/source/theory/mhd.rst`

- [ ] **Step 1: Read MHD implementation**

Read `src/core/numerics/` for MHD-related solvers. Search for reconstruction, Riemann solver, and MHD time stepper implementations. Also read `pyphare/pyphare/pharein/simulation.py` for valid option strings.

- [ ] **Step 2: Write mhd.rst**

Create new file with:

- Title: "The MHD formalism"
- **The ideal MHD equations** — conservation form with `.. math::`:
  - Mass: ∂ρ/∂t + ∇·(ρv) = 0
  - Momentum: ∂(ρv)/∂t + ∇·(ρvv + P*I - BB/μ₀) = 0
  - Energy: ∂E_tot/∂t + ∇·((E_tot + P*)v - B(v·B)/μ₀) = 0
  - Induction: ∂B/∂t - ∇×(v×B) = 0
  where P* = P + B²/(2μ₀) is the total pressure
- **Finite volume discretization** — brief explanation: integrate conservation laws over cells, fluxes at cell interfaces computed from reconstructed states via Riemann solver.
- **Reconstruction methods** — subsection describing each:
  - Constant (1st order, piecewise constant, most diffusive, most robust)
  - Linear + MinMod limiter (2nd order, slope-limited to prevent oscillations)
  - WENO3 (3rd order weighted ENO, nonlinear weights for shock capturing)
  - WENOZ (improved WENO weights, less dissipative than WENO3)
  - MP5 (5th order monotonicity-preserving, highest accuracy on smooth solutions)
- **Riemann solvers** — subsection:
  - Rusanov (local Lax-Friedrichs): uses maximum wave speed, most diffusive but never fails
  - HLLD: resolves all 7 MHD wave families, much less diffusive, recommended for production
- **Time stepping** — subsection:
  - Euler (1st order, for debugging only)
  - TVDRK2 (2nd order TVD Runge-Kutta, good balance of accuracy and cost)
  - SSPRK4_5 (4th order, 5-stage strong stability preserving, best for high-order spatial discretization)
- **Dissipative terms** — resistivity η, viscosity ν, hyper-resistivity. Same role as in hybrid PIC but applied to MHD equations.
- **Hall MHD** — brief note: when `hall=True`, the Hall term (J×B/ne) is added to the induction equation, capturing dispersive whistler waves.
- **Multi-model coupling** — PHARE can run MHD on coarse AMR levels and Hybrid PIC on fine levels. `max_mhd_level` sets the highest level using MHD. Above that level, the solver switches to PIC. This enables large-scale MHD context with kinetic physics at the regions of interest.

- [ ] **Step 3: Commit**

```bash
git add doc/source/theory/mhd.rst
git commit -m "docs: add MHD theory - conservation laws, reconstruction, Riemann solvers, time stepping, multi-model coupling"
```

---

## Task 9: Tutorial 1 — alfven_wave_1d.rst

**Files:**
- Create: `doc/source/tutorials/alfven_wave_1d.rst`

- [ ] **Step 1: Read source test for physics and parameters**

Read `tests/functional/alfven_wave/alfven_wave1d.py` to extract the exact configuration, initial conditions, and analysis approach.

- [ ] **Step 2: Write alfven_wave_1d.rst**

Create the file with:

- Title: "Tutorial: 1D Alfven Wave"
- **Introduction** — what this tutorial covers (first PHARE simulation), what an Alfven wave is (transverse MHD wave propagating along the magnetic field at speed VA = B₀/√(μ₀ n m)), what we expect to see (By oscillation propagating at VA).
- **The physics** — 2-3 paragraphs explaining the setup: uniform background B₀ along x, uniform density, small perturbation in By and Viy. The perturbation propagates as a circularly polarized Alfven wave.
- **Setting up the simulation** — full commented Python script:
  ```python
  import numpy as np
  import pyphare.pharein as ph

  # Physical parameters
  B0 = 1.0       # background magnetic field
  n0 = 1.0       # background density
  beta = 1.0     # plasma beta
  ...

  ph.Simulation(
      time_step=0.001,
      time_step_nbr=1000,
      cells=400,
      dl=0.3,
      ...
  )

  ph.MaxwellianFluidModel(
      bx=lambda x: B0,
      by=lambda x: ampl * np.cos(2 * np.pi * x / L),
      bz=lambda x: 0.0,
      protons={...}
  )

  ph.ElectronModel(closure="isothermal", Te=...)
  ph.ElectromagDiagnostics(quantity="B", write_timestamps=...)
  ```
  Walk through each parameter explaining *why* it's set that way.
- **Running the simulation** — show the execution command: `python3 -Ou alfven_wave_1d.py`
- **Analyzing results** — complete post-processing script:
  ```python
  from pyphare.pharesee.run import Run
  import matplotlib.pyplot as plt

  run = Run("./diags")
  for t in [0.0, 0.25, 0.5]:
      B = run.GetB(t, merged=True, interp="nearest")
      # extract and plot By
  ```
- **What to check** — the wave should propagate at the Alfven speed. Show how to measure the phase speed from the spatial shift of the By profile over time.
- Cross-references: :doc:`../usage/simulation_inputs`, :doc:`../pharesee/get_data`, :doc:`../theory/hybridpic`

- [ ] **Step 3: Commit**

```bash
git add doc/source/tutorials/alfven_wave_1d.rst
git commit -m "docs: add tutorial 1 - 1D Alfven wave (beginner)"
```

---

## Task 10: Tutorial 2 — harris_2d.rst

**Files:**
- Create: `doc/source/tutorials/harris_2d.rst`

- [ ] **Step 1: Read harris_2d.py source**

Read `tests/functional/harris/harris_2d.py` for configuration and physics.

- [ ] **Step 2: Write harris_2d.rst**

- Title: "Tutorial: 2D Harris Sheet Reconnection"
- **Introduction** — magnetic reconnection is a fundamental plasma process. The Harris current sheet is the standard test case.
- **The physics** — Harris equilibrium: Bx = B0 tanh(y/L), density = n0/cosh²(y/L) + n_bg. Explain why reconnection happens (tearing instability).
- **Setting up the simulation** — full 2D config:
  - `cells=(nx, ny)`, `dl=(dx, dy)`, `boundary_types=("periodic", "periodic")`
  - 2D initialization functions `f(x, y)` with the Harris profile
  - Static refinement boxes around y=Ly/2 (the current sheet)
  - Both EM and fluid diagnostics
- **Running** — `mpirun -np 4 python3 -Ou harris_2d.py` (2D simulations benefit from MPI)
- **Analyzing results** — 2D colormap of Jz at several times, showing reconnection developing. Show matplotlib `pcolormesh` example. Then show reconnection rate analysis:
  ```python
  Az, (xn, yn) = run.GetMagneticFlux(time)
  times_centered, rates, flux_at_xpoint, xpoint_traj = run.GetReconnectionRate(times)
  ```
- Cross-references: :doc:`../theory/amr`, :doc:`../pharesee/reconnection_analysis`

- [ ] **Step 3: Commit**

```bash
git add doc/source/tutorials/harris_2d.rst
git commit -m "docs: add tutorial 2 - 2D Harris sheet reconnection (intermediate)"
```

---

## Task 11: Tutorial 3 — multi_population_1d.rst

**Files:**
- Create: `doc/source/tutorials/multi_population_1d.rst`

- [ ] **Step 1: Read ion_ion_beam source**

Read `tests/functional/ionIonBeam/ion_ion_beam1d.py`.

- [ ] **Step 2: Write multi_population_1d.rst**

- Title: "Tutorial: Multi-Population Ion-Ion Beam"
- **Introduction** — counter-streaming ion beams drive electromagnetic instability. Demonstrates multi-population setup.
- **The physics** — two proton populations ("main" and "beam") with different bulk velocities. The relative drift drives instability and wave growth.
- **Setting up** — full config showing MaxwellianFluidModel with two populations:
  ```python
  ph.MaxwellianFluidModel(
      bx=..., by=..., bz=...,
      main={"charge": 1, "density": ..., "vbulkx": v_main, ...},
      beam={"charge": 1, "density": ..., "vbulkx": v_beam, ...},
  )
  ```
  Emphasize that population names are arbitrary keywords.
- **Per-population diagnostics** — show how to set up separate fluid and particle diagnostics for each population using `population_name`.
- **Analysis** — compare density profiles of each population, plot velocity distributions showing beam structure.
- Cross-reference: :doc:`../usage/simulation_inputs`

- [ ] **Step 3: Commit**

```bash
git add doc/source/tutorials/multi_population_1d.rst
git commit -m "docs: add tutorial 3 - multi-population ion-ion beam (intermediate)"
```

---

## Task 12: Tutorial 4 — mhd_shock_1d.rst

**Files:**
- Create: `doc/source/tutorials/mhd_shock_1d.rst`

- [ ] **Step 1: Read mhd_shock source**

Read `tests/functional/mhd_shock/mhd_shock.py`.

- [ ] **Step 2: Write mhd_shock_1d.rst**

- Title: "Tutorial: 1D MHD Shock Tube"
- **Introduction** — first MHD simulation. The Riemann problem is the fundamental test for any MHD solver.
- **The physics** — two constant states separated by a discontinuity. Evolution produces a rarefaction wave, contact discontinuity, and shock wave.
- **Setting up** — full MHD config:
  ```python
  ph.Simulation(
      model_options=["MHDModel"],
      reconstruction="WENOZ",
      riemann="Rusanov",
      mhd_timestepper="TVDRK2",
      ...
  )
  ph.MHDModel(
      density=lambda x: np.where(x < L/2, rho_L, rho_R),
      vx=lambda x: np.where(x < L/2, vx_L, vx_R),
      ...
  )
  ```
  Explain each MHD-specific parameter.
- **Analysis** — plot density, velocity, pressure profiles at final time. Compare left/right states.
- Cross-reference: :doc:`../theory/mhd`, :doc:`../usage/simulation_inputs`

- [ ] **Step 3: Commit**

```bash
git add doc/source/tutorials/mhd_shock_1d.rst
git commit -m "docs: add tutorial 4 - 1D MHD shock tube (MHD beginner)"
```

---

## Task 13: Tutorial 5 — mhd_orszag_tang_2d.rst

**Files:**
- Create: `doc/source/tutorials/mhd_orszag_tang_2d.rst`

- [ ] **Step 1: Read orszag_tang source**

Read `tests/functional/mhd_orszagtang/orszag_tang.py`.

- [ ] **Step 2: Write mhd_orszag_tang_2d.rst**

- Title: "Tutorial: 2D Orszag-Tang Vortex"
- **Introduction** — classic 2D MHD test. Smooth initial conditions develop into complex shock interactions.
- **The physics** — sinusoidal velocity and magnetic field initial conditions on a 2D periodic domain. The nonlinear evolution creates interacting shocks, current sheets, and vortices.
- **Setting up** — full 2D MHD config with SSPRK4_5 time stepping, WENOZ reconstruction, viscosity parameter.
- **Running** — note that 2D MHD at 256×256 is moderately expensive, recommend MPI.
- **Analysis** — 2D colormaps of density and pressure at several times showing shock formation. Qualitative comparison with published reference results.
- Cross-reference: :doc:`../theory/mhd`

- [ ] **Step 3: Commit**

```bash
git add doc/source/tutorials/mhd_orszag_tang_2d.rst
git commit -m "docs: add tutorial 5 - 2D Orszag-Tang vortex (MHD intermediate)"
```

---

## Task 14: Tutorial 6 — amr_tagging.rst

**Files:**
- Create: `doc/source/tutorials/amr_tagging.rst`

- [ ] **Step 1: Read tagging test source**

Read `tests/functional/translation/translat1d.py` for the tagged refinement configuration.

- [ ] **Step 2: Write amr_tagging.rst**

- Title: "Tutorial: Adaptive Mesh Refinement with Tagging"
- **Introduction** — static boxes require knowing where to refine a priori. Adaptive tagging lets PHARE decide at runtime based on the solution.
- **Static vs adaptive** — brief comparison. Static: `refinement="boxes"` with `refinement_boxes` dict. Adaptive: `refinement="tagging"` with `tagging_threshold`.
- **Setting up** — config with tagging enabled:
  ```python
  ph.Simulation(
      refinement="tagging",
      max_nbr_levels=3,
      tagging_threshold=0.5,
      nesting_buffer=1,
      ...
  )
  ```
  Explain each AMR parameter and how it affects refinement.
- **Monitoring refinement** — use MetaDiagnostics to output refinement tags:
  ```python
  ph.MetaDiagnostics(quantity="tags", write_timestamps=timestamps)
  ```
- **Analysis** — show how to load and visualize AMR levels: overlay patch boundaries on field data. Show how the refined regions track the physics (e.g., follow a propagating discontinuity).
- **Performance tips** — patch size trade-offs (too small → overhead, too large → waste), clustering algorithm choice.
- Cross-reference: :doc:`../theory/amr`, :doc:`../usage/simulation_inputs`

- [ ] **Step 3: Commit**

```bash
git add doc/source/tutorials/amr_tagging.rst
git commit -m "docs: add tutorial 6 - adaptive mesh refinement with tagging (advanced)"
```

---

## Task 15: Rewrite pharesee/get_data.rst

**Files:**
- Modify: `doc/source/pharesee/get_data.rst`

- [ ] **Step 1: Read Run class source for accurate API**

Read `pyphare/pyphare/pharesee/run/run.py` — check method signatures, return types, and keyword arguments.

- [ ] **Step 2: Rewrite get_data.rst**

Replace with comprehensive guide:

- Title: "Loading and accessing data"
- **Creating a Run object**:
  ```python
  from pyphare.pharesee.run import Run
  run = Run("/path/to/output_dir")
  ```
- **Querying available times**: `run.times("EM_B")` (key = HDF5 filename stem), `run.all_times()`
- **Getter methods** — organized tables by category:
  - Electromagnetic: `GetB(t)`, `GetE(t)`, `GetJ(t)`, `GetDivB(t)`
  - Ion fluid (total): `GetNi(t)`, `GetMassDensity(t)`, `GetVi(t)`, `GetPi(t)`, `GetPe(t)`
  - Ion fluid (per-population): `GetN(t, "protons")`, `GetFlux(t, "protons")`, `GetPressure(t, "protons")`
  - MHD: `GetMHDrho(t)`, `GetMHDV(t)`, `GetMHDP(t)`, `GetMHDrhoV(t)`, `GetMHDEtot(t)`
  - Particles: `GetParticles(t, "protons")`
  - Metadata: `GetTags(t)`, `GetRanks(t)`, `GetParticleCount(t)`
- **Common keyword arguments** — explain each:
  - `merged=True/False`: when True, interpolates all AMR levels onto finest grid
  - `interp="nearest"` or `"linear"`: interpolation method for merging
  - `all_primal=True`: moves staggered data to primal grid nodes for easier plotting
  - `selection_box`: restrict query to a spatial subregion
- **The PatchHierarchy** — what the return value looks like:
  - `hier.levels(time)` → dict of PatchLevel
  - Iterating: `for patch in level.patches`
  - Accessing data: `patch.patch_datas["By"].dataset`, `.x` coordinates
  - `hier.quantities()` — list available field names
- **Working with merged data** — show complete example of extracting 1D field values
- **Direct file loading**: `hierarchy_from("file.h5", time=t)` for advanced use
- Keep the `.. autoclass::` directive for Run at the bottom

- [ ] **Step 3: Commit**

```bash
git add doc/source/pharesee/get_data.rst
git commit -m "docs: comprehensive rewrite of pharesee data loading guide"
```

---

## Task 16: Fill pharesee/plotting_fields.rst

**Files:**
- Modify: `doc/source/pharesee/plotting_fields.rst`

- [ ] **Step 1: Read plotting utilities**

Read `pyphare/pyphare/pharesee/plotting.py` and check for `plot1d`, `plot2d` method signatures in the hierarchy classes.

- [ ] **Step 2: Write plotting_fields.rst**

Replace stub with:

- Title: "Plotting fields"
- **1D field plotting** — complete matplotlib example:
  ```python
  from pyphare.pharesee.run import Run
  import matplotlib.pyplot as plt

  run = Run("diags")
  B = run.GetB(0.5, merged=True, interp="nearest")

  # B is a PatchHierarchy. For 1D merged data:
  for patch in B.level(0).patches:
      by = patch.patch_datas["By"]
      plt.plot(by.x, by.dataset)
  plt.xlabel("x")
  plt.ylabel("By")
  plt.show()
  ```
- **Plotting time evolution** — multiple times on same axes with legend
- **2D field plotting** — `pcolormesh` example for 2D data
- **Contour plots** — for magnetic flux (vector potential)
- **Convenience methods** — `hier.plot1d()` and `hier.plot2d()` with their parameters
- **AMR patch boundaries** — how to overlay patch boundaries on 2D plots
- **Multi-panel figures** — side-by-side Bx, By, Bz using `plt.subplots`

- [ ] **Step 3: Commit**

```bash
git add doc/source/pharesee/plotting_fields.rst
git commit -m "docs: fill field plotting guide with 1D, 2D, and multi-panel examples"
```

---

## Task 17: Fill pharesee/plotting_distributions.rst

**Files:**
- Modify: `doc/source/pharesee/plotting_distributions.rst`

- [ ] **Step 1: Read particle data structures**

Read how particle data is stored in the hierarchy — check `pyphare/pyphare/pharesee/hierarchy/` for particle patch data classes.

- [ ] **Step 2: Write plotting_distributions.rst**

Replace stub with:

- Title: "Plotting particle distributions"
- **Loading particle data**:
  ```python
  particles = run.GetParticles(time, "protons")
  ```
- **Particle attributes** — explain the data structure:
  - `.iCell` + `.delta` give position (cell index + offset within cell)
  - `.v` gives velocity (3-component, always 3D even in 1D sims)
  - `.weight` gives statistical weight of each macro-particle
  - `.charge` gives charge
- **Phase space plots** — x vs vx scatter plot (weighted):
  ```python
  # compute physical position from iCell and delta
  x = (particles.iCell[:, 0] + particles.delta[:, 0]) * dx
  vx = particles.v[:, 0]
  plt.scatter(x, vx, s=0.1, alpha=0.1)
  ```
- **1D velocity distributions** — weighted histogram of vx:
  ```python
  plt.hist(particles.v[:, 0], bins=100, weights=particles.weight)
  ```
- **2D velocity space** — vx-vy 2D histogram using `plt.hist2d` with weights
- **Convenience method** — `hier.dist_plot()` if available
- **Tips** — downsample for large particle counts; always use weights; consider density-weighted histograms for physical distributions

- [ ] **Step 3: Commit**

```bash
git add doc/source/pharesee/plotting_distributions.rst
git commit -m "docs: fill particle distribution plotting guide with phase space and velocity histograms"
```

---

## Task 18: Create pharesee/reconnection_analysis.rst

**Files:**
- Create: `doc/source/pharesee/reconnection_analysis.rst`

- [ ] **Step 1: Read reconnection analysis methods**

Read `pyphare/pyphare/pharesee/run/run.py` — look at `GetMagneticFlux`, `FindPrimaryXPoint`, `GetReconnectionRate` method signatures and return values.

- [ ] **Step 2: Write reconnection_analysis.rst**

- Title: "Reconnection analysis"
- **Introduction** — specialized tools for analyzing 2D magnetic reconnection simulations. Only applicable to 2D runs with magnetic field diagnostics.
- **Computing the vector potential**:
  ```python
  Az, (xn, yn) = run.GetMagneticFlux(time, interp="nearest")
  ```
  Explain: Az is the z-component of the vector potential, computed by integrating B. xn, yn are the coordinate grids.
- **Finding X-points**:
  ```python
  x_xpoint, y_xpoint, idx = run.FindPrimaryXPoint(Az, xn, yn)
  ```
  Returns a 3-tuple: x-coordinate, y-coordinate, and array index of the primary X-point (the saddle point in Az closest to the center of the domain).
- **Reconnection rate**:
  ```python
  times = [0.0, 10.0, 20.0, 30.0]
  times_centered, rates, flux_at_xpoint, xpoint_trajectory = run.GetReconnectionRate(times)
  ```
  Returns a 4-tuple: centered times, reconnection rates (dAz/dt at X-point), magnetic flux values at the X-point, and X-point position trajectory.
- **Worked example** — complete script from a Harris sheet simulation: load data → compute Az → find X-point → plot contours of Az with X-point marked → plot reconnection rate vs time
- Cross-reference: :doc:`../tutorials/harris_2d`

- [ ] **Step 3: Commit**

```bash
git add doc/source/pharesee/reconnection_analysis.rst
git commit -m "docs: add reconnection analysis guide - vector potential, X-points, reconnection rate"
```

---

## Task 19: Final review and Sphinx build test

- [ ] **Step 1: Verify all RST files exist**

```bash
ls -la doc/source/theory/spatial_discretization.rst \
      doc/source/theory/temporal_discretization.rst \
      doc/source/theory/amr.rst \
      doc/source/theory/mhd.rst \
      doc/source/tutorials/alfven_wave_1d.rst \
      doc/source/tutorials/harris_2d.rst \
      doc/source/tutorials/multi_population_1d.rst \
      doc/source/tutorials/mhd_shock_1d.rst \
      doc/source/tutorials/mhd_orszag_tang_2d.rst \
      doc/source/tutorials/amr_tagging.rst \
      doc/source/pharesee/reconnection_analysis.rst
```

- [ ] **Step 2: Check cross-references**

Grep all RST files for `:doc:` and `:ref:` directives. Verify each target file exists and the path is correct relative to the referencing file.

```bash
grep -rn ':doc:' doc/source/ | grep -v '_build'
```

- [ ] **Step 3: Check that usage/examples.rst is gone**

```bash
test ! -f doc/source/usage/examples.rst && echo "OK: examples.rst removed"
```

- [ ] **Step 4: Attempt Sphinx build (if dependencies available)**

```bash
cd doc && make html 2>&1 | tail -20
```

If Sphinx is not installed, skip this step — it's a nice-to-have.

- [ ] **Step 5: Final commit if any fixes needed**

```bash
git add doc/source/
git commit -m "docs: fix any remaining cross-references and build warnings"
```

Only create this commit if there are actual changes to commit.
