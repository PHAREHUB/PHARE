# PHARE End-User Documentation Design

**Date**: 2026-03-17
**Status**: Approved
**Format**: Sphinx RST, filling existing `doc/source/` structure
**Audience**: All levels — from grad students new to Hybrid PIC to experienced PHARE users

---

## Overview

Comprehensive rewrite and expansion of PHARE's Sphinx documentation. The existing structure has good content in `theory/hybridpic.rst`, `theory/pic.rst`, and `usage/simulation_inputs.rst`, but most pages are empty stubs. This spec fills every gap and adds new sections for MHD (equal treatment with Hybrid PIC), progressive tutorials, and data analysis.

**21 pages total**: 6 existing kept/expanded, 4 empty pages filled, 11 new pages.

---

## 1. Theory Section

### 1.1 Keep as-is
- `theory/hybridpic.rst` — Vlasov equation, Ohm's law, quasineutrality, closures, dissipative terms
- `theory/pic.rst` — Macro-particles, b-spline shape functions, interpolation orders

### 1.2 Fill: `theory/spatial_discretization.rst`
**Content:**
- Yee lattice layout: primal vs dual centering for each field component
- How centering differs per dimension (1D/2D/3D)
- Finite difference stencils for curl operators (Faraday: curl E, Ampere: curl B)
- Gradient operator (for electron pressure term in Ohm's law)
- Ghost cells: purpose, width depends on interpolation order
- Diagrams showing field placement on the staggered grid in 1D and 2D

### 1.3 Fill: `theory/temporal_discretization.rst`
**Content:**
- The predictor-corrector scheme (SolverPPC algorithm):
  1. predictor1: Faraday(B→Bpred), Ampere(Bpred→J), Ohm(→Epred)
  2. average: Bavg=(B+Bpred)/2, Eavg=(E+Epred)/2
  3. moveIons: push domain particles with averaged fields
  4. predictor2: repeat with updated moments
  5. average: update averages
  6. moveIons: push all particles (domain + ghosts)
  7. corrector: Faraday(→B_final), Ampere(→J_final), Ohm(→E_final)
- Boris pusher algorithm: half-acceleration → rotation (t, s vectors) → half-acceleration → position update
- Stability constraints: particles must not cross more than one cell per timestep
- Subcycling: finer AMR levels take smaller timesteps proportional to the **square** of the refinement ratio (dt_fine = dt_coarse / ratio^2). With refinement ratio 2, this means 4 substeps per coarse step. This is because whistler waves impose dt ~ dx^2 scaling.

### 1.4 Fill: `theory/amr.rst`
**Content:**
- What AMR is: concentrate resolution where needed, save computation elsewhere
- Patch-based structured AMR (SAMRAI library)
- Hierarchy: levels, patches, cells. Level 0 = coarsest, level N = finest
- Refinement ratio: currently fixed at 2 in PHARE (not user-configurable). Cell size halves at each level.
- Two refinement strategies:
  - Static boxes: user specifies regions a priori
  - Adaptive tagging: B-field gradient criterion, `tagging_threshold`
- Ghost filling between levels:
  - Spatial interpolation from coarser to finer level
  - Time interpolation for subcycling (coarse level data at intermediate times)
- Particles at level boundaries:
  - `refined_particle_nbr`: how many fine particles replace one coarse particle
  - Level ghost particles: filled from coarser level via splitting
  - Patch ghost particles: filled from neighboring patches at same level
- Patch size constraints: `smallest_patch_size`, `largest_patch_size`, `nesting_buffer`
- Clustering algorithms: Berger-Rigoutsos vs tile-based (default: "tile")

### 1.5 New: `theory/mhd.rst`
**Content:**
- Ideal MHD equations in conservation form (mass, momentum, energy, induction)
- Finite volume discretization
- Reconstruction methods available in PHARE:
  - Constant (1st order, most diffusive)
  - Linear with limiters (2nd order): MinMod
  - WENO3 (3rd order weighted ENO)
  - WENOZ (improved WENO weights)
  - MP5 (5th order monotonicity-preserving)
- Riemann solvers:
  - Rusanov (most diffusive, most robust)
  - HLLD (less diffusive, handles all MHD waves)
- Time stepping schemes:
  - Euler (1st order)
  - TVDRK2 (2nd order TVD Runge-Kutta)
  - SSPRK4_5 (4th order strong stability preserving, 5 stages)
- Dissipation: resistivity (eta), viscosity (nu), hyper-resistivity
- Hall MHD extension
- Multi-model coupling: MHD on coarse levels, Hybrid PIC on fine levels (`max_mhd_level`)

---

## 2. Getting Started Section

### 2.1 Expand: `getting.rst`
**Content:**
- Prerequisites table:
  - C++20 compiler: gcc >= 10 or clang >= 13
  - CMake >= 3.20.1
  - MPI implementation (OpenMPI, MPICH, Intel MPI)
  - Parallel HDF5 (built with MPI support)
  - Python >= 3.8
- Clone command with `--recursive` for submodules
- Auto-downloaded dependencies (users do NOT need to install): SAMRAI, pybind11, HighFive, GoogleTest
- Python dependencies: `pip install -r requirements.txt`
- Virtual environment recommendation

### 2.2 Expand: `build.rst`
**Content:**
- Production build:
  ```
  mkdir build && cd build
  cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-O3 -march=native"
  cmake --build . -j$(nproc)
  ```
- Debug build with symbols
- CMake options table (all flags with defaults and descriptions):
  - `devMode`, `test`, `testMPI`, `PHARE_MPI_PROCS`
  - `asan`, `ubsan`, `coverage`
  - `bench`, `SAMRAI_ROOT`
  - `lowResourceTests`, `PHARE_EXEC_LEVEL_MIN/MAX`
- Setting PYTHONPATH (critical):
  ```
  export PYTHONPATH=/path/to/PHARE/pyphare:/path/to/build:$PYTHONPATH
  ```
- Verifying the build: smoke test script that imports the pybind module and runs a trivial simulation
- Common build issues and solutions:
  - HDF5 not parallel → check `h5cc -showconfig`
  - Wrong Python picked up → use `-DPython_EXECUTABLE=...`
  - MPI not found → set `MPI_HOME`
  - SAMRAI download fails → use `SAMRAI_ROOT` for pre-built

### 2.3 Expand: `usage/run_from_python.rst`
**Content (keep existing, fix typos, add):**
- Fix "minimmum" → "minimum"
- MPI execution: `mpirun -np 4 python3 -Ou my_script.py`
- Explain `-Ou` flag (unbuffered output)
- Simulator lifecycle: setup() → initialize() → run()/advance() → reset()
- Post-advance callback usage
- Log files (per-rank logs when `log_to_file=True`)
- Dry run mode: `dry_run=True` or `PHARE_DRY_RUN=1`

### 2.4 Stub: `usage/run_from_exe.rst`
- Note that Python is the only supported entry point for PHARE
- Redirect to `run_from_python.rst`

---

## 3. Usage Section (Simulation Configuration)

### 3.1 Rewrite: `usage/simulation_inputs.rst`
Complete reference for all configuration blocks. Fix existing typos. Organize by subsections:

**3.1.1 The Simulation block**
All parameters organized by category with types, defaults, and descriptions:
- Domain & mesh: `cells`, `dl`, `domain_size`, `boundary_types`
- Time: `time_step`, `time_step_nbr`, `final_time` (2-of-3 rule)
- Numerics: `interp_order`, `particle_pusher` (default "modified_boris"), `layout`, `refined_particle_nbr`
- AMR: `max_nbr_levels`, `refinement`, `refinement_boxes`, `tagging_threshold`, `smallest_patch_size`, `largest_patch_size`, `nesting_buffer`, `tag_buffer`, `clustering`
- Physics: `resistivity` (default 0.0), `hyper_resistivity` (default 0.0001 — note: non-zero by default), `hyper_mode`
- MHD: `model_options`, `gamma`, `eta`, `nu`, `hall`, `res`, `hyper_res`, `reconstruction`, `limiter`, `riemann`, `mhd_timestepper`, `max_mhd_level`
- Diagnostics: `diag_options` dict structure
- Restart: `restart_options` dict structure
- Misc: `strict`, `dry_run`, `description`, `write_reports`

**3.1.2 MaxwellianFluidModel**
- Magnetic field: `bx`, `by`, `bz` as callables
- Ion populations as keyword arguments with dict values:
  - `charge`, `mass`, `density`, `vbulkx/y/z`, `vthx/y/z`, `nbr_part_per_cell`
  - Optional: `init` (with `seed`), `density_cut_off`
- Function signatures per dimension: `f(x)`, `f(x,y)`, `f(x,y,z)`
- Functions must be periodic (validated at boundaries)
- Multi-population example
- Keep existing Harris sheet example (it's good)

**3.1.3 MHDModel**
- `density`, `vx` (default 1.0 — note asymmetry), `vy` (default 0.0), `vz` (default 0.0), `bx`, `by`, `bz`, `p` as spatial callables
- Must set `model_options=["MHDModel"]` in Simulation
- Guidance on choosing reconstruction/limiter/Riemann solver:
  - Start with WENOZ + Rusanov for robustness
  - Use HLLD for less diffusive results once stable
  - MP5 for highest accuracy on smooth problems

**3.1.4 ElectronModel**
- `closure="isothermal"` (only option currently)
- `Te`: electron temperature in simulation units (energy / density)
- Physical meaning: Pe = n * Te everywhere

**3.1.5 Diagnostics**
- Common parameters: `write_timestamps`, `flush_every`, `quantity`
  - Note: `elapsed_timestamps` and `attributes` exist in some code paths but may not pass validation in all cases — document `write_timestamps` as the primary interface
- ElectromagDiagnostics: `quantity` = "E" or "B"
- FluidDiagnostics:
  - Total ion (no `population_name` needed): "charge_density", "mass_density", "bulkVelocity", "momentum_tensor"
  - Per-population (requires `population_name`): "flux", "momentum_tensor", "density"
  - Special: "pressure_tensor" — convenience shorthand that auto-creates sub-diagnostics (density + flux + momentum_tensor for the population). Document that this is a high-level wrapper, not a direct quantity.
- ParticleDiagnostics: `quantity` = "domain", "levelGhost", or "space_box" (with `extent`)
  - Requires `population_name`
  - Warning: particle output is much larger than field output
- MHDDiagnostics: `quantity` = "rho", "V", "P", "rhoV", "Etot"
- MetaDiagnostics: `quantity` = "tags"
- InfoDiagnostics: `quantity` = "particle_count"
- Timestamp generation patterns (uniform spacing, every-N-steps)
- Diagnostic output formats: "phareh5" (HDF5, default) and "pharevtkhdf" (VTK+HDF5, useful for 3D visualization)

**3.1.6 LoadBalancer**
- `active`, `mode` ("nppc" / "homogeneous"), `tol`, `auto`, `every`, `on_init`
- Expert parameters (when `auto=True`): `next_rebalance`, `next_rebalance_backoff_multiplier`, `max_next_rebalance`
- When to use: large particle count imbalance across MPI ranks
- nppc mode: balance by particle count (recommended for PIC)
- homogeneous mode: balance by cell count

---

## 4. Tutorials Section (New)

All tutorials are self-contained scripts with physics explanation, parameter walkthrough, and post-processing.

### 4.1 New: `tutorials/alfven_wave_1d.rst`
**Level**: Beginner — first PHARE simulation
**Physics**: Circularly polarized Alfven wave propagating in a uniform magnetized plasma
**What the user learns**:
- Simulation block (1D, periodic, time parameters)
- MaxwellianFluidModel with one population
- ElectronModel
- ElectromagDiagnostics
- Running with Simulator
- Loading results with Run.GetB()
- Plotting By(x) at multiple times
- Measuring phase speed and comparing with Alfven speed VA = B/sqrt(mu0 * n * m)

**Source inspiration**: `tests/functional/alfven_wave/alfven_wave1d.py`

### 4.2 New: `tutorials/harris_2d.rst`
**Level**: Intermediate — 2D simulation with AMR
**Physics**: Harris current sheet equilibrium (Bx reversal, density enhancement at sheet), tearing instability, magnetic reconnection
**What the user learns**:
- 2D domain setup: `cells=(nx, ny)`, `dl=(dx, dy)`
- 2D initialization functions: `f(x, y)`
- Static refinement boxes around the current sheet
- FluidDiagnostics and ParticleDiagnostics
- Loading and plotting 2D fields (Jz colormap)
- Computing reconnection rate with `GetReconnectionRate`

**Source inspiration**: `tests/functional/harris/harris_2d.py`

### 4.3 New: `tutorials/multi_population_1d.rst`
**Level**: Intermediate — multiple ion species
**Physics**: Ion-ion beam instability — two proton populations counter-streaming, generating electromagnetic instability
**What the user learns**:
- Defining multiple populations in MaxwellianFluidModel
- Per-population diagnostics with `population_name`
- Comparing density and velocity of each population
- Phase space plots (x vs vx) for each population

**Source inspiration**: `tests/functional/ionIonBeam/ion_ion_beam1d.py`

### 4.4 New: `tutorials/mhd_shock_1d.rst`
**Level**: Beginner (MHD) — first MHD simulation
**Physics**: 1D Riemann problem (Sod-like shock tube) — discontinuous initial conditions evolve into rarefaction + contact + shock
**What the user learns**:
- MHDModel instead of MaxwellianFluidModel
- `model_options=["MHDModel"]` in Simulation
- Reconstruction, limiter, Riemann solver choices
- MHDDiagnostics (rho, V, P)
- Post-processing: density/velocity/pressure profiles
- Comparison with analytical solution (optional)

**Source inspiration**: `tests/functional/mhd_shock/mhd_shock.py`

### 4.5 New: `tutorials/mhd_orszag_tang_2d.rst`
**Level**: Intermediate (MHD) — complex 2D MHD flow
**Physics**: Orszag-Tang vortex — smooth initial conditions develop interacting shocks, showcases code's ability to handle complex MHD dynamics
**What the user learns**:
- 2D MHD setup with larger grids
- Viscosity parameter (`nu`)
- Higher-order time stepping (SSPRK4_5)
- 2D colormap plots of density and pressure at several times
- Qualitative comparison with published results

**Source inspiration**: `tests/functional/mhd_orszagtang/orszag_tang.py`

### 4.6 New: `tutorials/amr_tagging.rst`
**Level**: Advanced — dynamic adaptive mesh refinement
**Physics**: Tangential discontinuity or Harris sheet with adaptive tagging — the grid refines dynamically as structures develop
**What the user learns**:
- `refinement="tagging"` vs `refinement="boxes"`
- `tagging_threshold`, `nesting_buffer`, `max_nbr_levels`
- MetaDiagnostics for refinement tags
- Visualizing AMR levels overlaid on field data
- How patch boundaries and levels evolve in time
- Performance considerations: patch size, clustering algorithm

**Source inspiration**: `tests/functional/translation/translat1d.py` (td1dtagged)

---

## 5. Data Analysis Section (pharesee)

### 5.1 Rewrite: `pharesee/get_data.rst`
**Content:**
- Creating a `Run` object: `Run("/path/to/output")`
- Querying available data: `run.times("EM_B")` (key matches HDF5 filename stem), `run.all_times()`
- Getter methods organized by category:
  - Electromagnetic: `GetB`, `GetE`, `GetJ`, `GetDivB`
  - Ion fluid (total): `GetNi`, `GetMassDensity`, `GetVi`, `GetPi`, `GetPe`
  - Ion fluid (per-population): `GetN(t, pop)`, `GetFlux(t, pop)`, `GetPressure(t, pop)`
  - MHD: `GetMHDrho`, `GetMHDV`, `GetMHDP`, `GetMHDrhoV`, `GetMHDEtot`
  - Particles: `GetParticles(t, pop)`
  - Metadata: `GetTags`, `GetRanks`, `GetParticleCount`
- Common keyword arguments:
  - `merged=True`: interpolate all levels onto finest grid
  - `interp="nearest"` or `"linear"`: interpolation method
  - `all_primal=True`: move staggered data to primal nodes
  - `selection_box`: restrict to spatial subregion
- The PatchHierarchy object:
  - `hier.levels(time)` → dict of PatchLevel
  - `hier.level(lvl_num, time)` → PatchLevel
  - `hier.quantities()` → available field names
  - Iterating patches: `for patch in level.patches`
  - Accessing data: `patch.patch_datas["Bx"].dataset`, `.x`, `.y`
- Direct loading: `hierarchy_from("file.h5", time=t)` (auto-dispatches between h5 and vtkhdf formats). Format-specific: `hierarchy_fromh5()` for HDF5 files.

### 5.2 Fill: `pharesee/plotting_fields.rst`
**Content:**
- 1D field plotting:
  - Extract coordinates and values from a merged hierarchy
  - matplotlib example: `plt.plot(x, by)`
  - Multiple times on same axes for evolution
- 2D field plotting:
  - Colormaps with `plt.pcolormesh` or `plt.imshow`
  - Contour plots for magnetic flux
  - Colorbar labeling with physical units
- Convenience methods: `hier.plot1d()`, `hier.plot2d()`
- Overlaying AMR patch boundaries on 2D plots
- Multi-panel figures (e.g., Bx, By, Bz side by side)

### 5.3 Fill: `pharesee/plotting_distributions.rst`
**Content:**
- Loading particle data: `particles = run.GetParticles(time, "protons")`
- Particle attributes: position from `iCell + delta`, velocity `.v`, `.weight`, `.charge`
- Phase space plots: x vs vx scatter (weighted by particle weight)
- 1D velocity distributions: histogram of vx (weighted)
- 2D velocity space: vx-vy 2D histogram
- Convenience method: `hier.dist_plot()`
- Tips: downsample for large particle counts, use density-weighted histograms

### 5.4 New: `pharesee/reconnection_analysis.rst`
**Content:**
- Computing vector potential: `Az, (xn, yn) = run.GetMagneticFlux(time)` (note: returns Az and a tuple of coordinates)
- Finding X-points: `run.FindPrimaryXPoint(Az, xn, yn)`
- Reconnection rate time series: `run.GetReconnectionRate(times)`
- Worked example from a 2D Harris sheet simulation
- Plotting: contour plot of Az with X-point marked, reconnection rate vs time

---

## 6. Updated `index.rst`

```rst
.. toctree::
   :caption: THEORY
   :maxdepth: 1
   :hidden:

   theory/hybridpic
   theory/pic
   theory/spatial_discretization
   theory/temporal_discretization
   theory/amr
   theory/mhd


.. toctree::
   :caption: BUILDING
   :maxdepth: 1
   :hidden:

   getting
   build


.. toctree::
   :caption: USAGE
   :maxdepth: 1
   :hidden:

   usage/simulation_inputs
   usage/run_from_python
   usage/run_from_exe


.. toctree::
   :caption: TUTORIALS
   :maxdepth: 1
   :hidden:

   tutorials/alfven_wave_1d
   tutorials/harris_2d
   tutorials/multi_population_1d
   tutorials/mhd_shock_1d
   tutorials/mhd_orszag_tang_2d
   tutorials/amr_tagging


.. toctree::
   :caption: ANALYSIS
   :maxdepth: 1
   :hidden:

   pharesee/get_data
   pharesee/plotting_fields
   pharesee/plotting_distributions
   pharesee/reconnection_analysis
```

---

## File Inventory

| File | Action | Priority |
|------|--------|----------|
| `theory/hybridpic.rst` | Keep | - |
| `theory/pic.rst` | Keep | - |
| `theory/spatial_discretization.rst` | Fill | High |
| `theory/temporal_discretization.rst` | Fill | High |
| `theory/amr.rst` | Fill | High |
| `theory/mhd.rst` | New | High |
| `getting.rst` | Expand | High |
| `build.rst` | Expand | High |
| `usage/simulation_inputs.rst` | Rewrite | High |
| `usage/run_from_python.rst` | Expand | Medium |
| `usage/run_from_exe.rst` | Stub | Low |
| `usage/examples.rst` | Remove (replaced by tutorials/) | - |
| `tutorials/alfven_wave_1d.rst` | New | High |
| `tutorials/harris_2d.rst` | New | High |
| `tutorials/multi_population_1d.rst` | New | Medium |
| `tutorials/mhd_shock_1d.rst` | New | High |
| `tutorials/mhd_orszag_tang_2d.rst` | New | Medium |
| `tutorials/amr_tagging.rst` | New | Medium |
| `pharesee/get_data.rst` | Rewrite | High |
| `pharesee/plotting_fields.rst` | Fill | Medium |
| `pharesee/plotting_distributions.rst` | Fill | Medium |
| `pharesee/reconnection_analysis.rst` | New | Low |
| `index.rst` | Update toctree | High |

---

## Writing Guidelines

- Write for all levels: brief intuitive explanation first, then precise details
- Every code block must be a complete, runnable snippet (no `# ...` placeholders)
- Use consistent import style: `import pyphare.pharein as ph`
- Parameter tables with columns: Name, Type, Default, Description
- Cross-reference between pages using `:doc:` and `:ref:`
- Fix all existing typos in kept/expanded pages
- No top-level imports of h5py, mpi4py, scipy — import at function scope in examples
- RST conventions: use `.. code-block:: python` for code, `.. math::` for equations, `.. note::` for tips, `.. warning::` for gotchas
