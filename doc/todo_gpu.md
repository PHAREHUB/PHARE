# GPU TODO: shortcuts to revisit

Places where the GPU (`GPU_UNIFIED`) path currently takes a shortcut that works but is not the
intended final design. Each entry says where, what the shortcut is, and what a proper fix looks
like.

## Particles

### Splitter partitions refined particles on the host
- `src/amr/data/particles/refine/particles_data_split.hpp` (`ParticlesRefining::forBoxes`)
- Always uses `ParticleArrayPartitioner<AllocatorMode::CPU, ...>`, whatever the simulation
  allocator. `forBoxes` runs on the host, per particle, over a small on-stack `AoSArray`, so the
  thrust GPU partitioner would launch a kernel per particle on host stack memory.
- Proper fix: do the whole split + partition for a patch on the device in one pass.

### AoSPCTS particle selection runs on the host
- `src/core/data/particles/selecting/detail/aos_selecting.hpp`
  (`ParticlesSelector<AoSPCTS, GPU_UNIFIED>::select/count`)
- Forward to the CPU implementations, reading `GPU_UNIFIED` memory from the host.
- Proper fix: device selection (count + scan + copy) per tile.

### AoS into AoSPC appending runs on the host
- `src/core/data/particles/appending/detail/aos_appending.hpp`
  (`ParticlesAppender<AoS, GPU_UNIFIED, AoSPC, GPU_UNIFIED>`)
- Forwards to the CPU `emplace_back` loop over `GPU_UNIFIED` memory.
- Proper fix: device append (per cell counts + scan + copy).

### AoSCMTS is CPU only
- `res/sim/gpu.txt`, `tests/simulator/layouts/__init__.py` (`allocators()`)
- The per-tile cell map is host side and must be re-registered on the vector after a move
  (`multi_boris.hpp`, "a cellmap can't be updated through its view").
- Proper fix: device-side cell map maintenance, or drop AoSCMTS for GPU entirely.

### AoSPCTS has no fused push + deposit GPU kernel
- `src/core/numerics/pusher/boris/detail/multi_boris.hpp` (`has_gpu_copy`)
- Only flat tiles (AoSTS, AoSMapped) have `on_tiles_push_deposit`, AoSPCTS goes through
  `move_rest` with a separate deposit.

### Flat tile GPU kernels are never exercised
- `src/core/numerics/pusher/boris/detail/multi_boris.hpp` (`on_tiles_push_deposit`),
  `src/core/numerics/ion_updater/ion_updater.hpp` (`mkn_xyz::ParallelIonUpdater`, flat deposit)
- No simulation, test or bench builds AoSTS/AoSMapped with `GPU_UNIFIED`, so the fused push +
  deposit kernel and the flat tile deposit are not even compiled for the device.
- Proper fix: AoSTS `GPU_UNIFIED` permutations in `test_multi_updater.cpp`.

### CFL check on the GPU: done except the fused flat kernel
- `src/core/numerics/pusher/boris/basics.hpp` (`MoveError`, `advance(p, halfdt, error)`),
  `multi_boris.hpp` (`MultiBoris::move_error()/check_move_error()`, functor member `move_error`),
  `ion_updater.hpp` (every updater calls `check_move_error()` after its moves / stream joins)
- With an error record `advance` records the first |delta| > 2 (any allocator, `compare_and_swap`)
  and returns the particle's own cell; the host throws a `DictionaryException` after the step.
  Without a record the CPU path still throws `MoveTwoCellException`.
- Left: `on_tiles_push_deposit` (static, flat tiles, HW only) does not get the record yet.
- Test: `tests/core/numerics/pusher/test_boris_move_error.cpp` (4/4 on ROCm, incl. a kernel).
  Not yet run under the mkn CPU emulation (`mkn.sh`), nor the CMake build.

### CPU emulation runs no shared memory kernels
- `src/core/numerics/pusher/boris/detail/multi_boris.hpp` (`on_tiles_push_deposit`),
  `src/core/numerics/interpolator/interpolating.hpp` (`on_tiles` and friends)
- With mkn's CPU fallthrough (`MKN_GPU_CPU=1`, `mkn.sh`) the AoSPCTS `GPU_UNIFIED` kernels and the
  field solvers run on the pseudo GPU (`test_multi_updater`, `test_field_solvers`), but the kernels
  using `extern __shared__` stay behind `PHARE_HAVE_MKN_GPU_HW`: there is no dynamic shared memory
  on the pseudo GPU.
- Proper fix: a per launch shared memory buffer in mkn.gpu's CPU backend, sized from the launch's
  `ds`, and kernels taking their shared memory from it instead of `extern __shared__`.

## Tiling

### 3D GPU tile size is a fixed default
- `src/core/data/tiles/tiling_options.hpp` (`TilingOptions::defaults_for`)
- `GPU_UNIFIED` in 3D defaults to `max_tile_size = 6` so a tile's ghost box of field nodes fits
  the GPU kernels' shared memory (`ds * 4 * 8` + warp buffer, checked against the device limit at
  launch, see `DeviceLauncher::check_shared_mem_`). Sized for interp order 1 only, interp >= 2 in
  3D does not fit.
- Proper fix: derive the max tile size from the ghost width (interp order) and the device's
  shared memory limit.

## Boundary conditions

### Tiled physical boundary conditions gather/scatter a whole patch
- `src/amr/data/field/refine/field_refine_patch_strategy.hpp` (`setPhysicalBoundaryConditions`)
- Tiled fields are reduced into a temporary untiled grid, the boundary condition is applied there,
  and the result is scattered back into every tile. Allocates per call and is host only.
- Proper fix: apply boundary conditions per tile on the tile's ghost region.

## Session state 2026-10-28 (uncommitted on tiling_gpu)

- `throw_exception` (`src/core/def.hpp`) now branches on `PHARE_DEVICE_PASS`: before, every HIP/CUDA
  build compiled host throws to `PHARE_ASSERT(false)`, silently disabling them.
- `IonUpdaterFor_t` (end of `ion_updater.hpp`) now picks `mkn_xyz::ParallelIonUpdater` for
  `GPU_UNIFIED` (not AoSMapped) when `PHARE_HAVE_MKN_GPU`; the host choice is
  `detail::HostIonUpdaterFor_t`; `solver_ppc.hpp` just uses `IonUpdaterFor_t`.
- Verified: ROCm `test_multi_updater` 42/42 and `test_boris_move_error` 4/4. To run: the emulation
  build (`mkn.sh`, needs `bin/core` rebuilt for CPU) for both tests, and a CMake build for
  `solver_ppc.hpp` / the new pusher test target.
- Seen once: a huge-dt `updatePopulations` (dt=1e3) on AoSPCTS CPU `parallel` died silently instead
  of throwing (test since removed); worth a look.

## Housekeeping

Not GPU shortcuts, but found reviewing #1416 and not to be forgotten.

### Dependencies unpinned and backends enabled by header presence
- `res/cmake/dep/*` (mkn.kul, mkn.gpu from master, Kokkos and KokkosTools from develop, magic_enum
  from HEAD), `src/core/def/kokkos_tools.hpp` and friends (`__has_include`)
- An installed Kokkos, RAJA, Umpire or mkn.gpu header turns that backend on whatever the CMake
  options say.
- Proper fix: pin each dependency to a tag/commit, enable backends from CMake options only.

### Unregistered tests
- 29 test sources are built by neither CMake nor `mkn.yaml`: the GPU ones (`test_*gpu*.cpp`,
  `test_raja_main.cpp`, `test_gpu_atomics*.cpp`, `test_gpu_map.cpp`, `test_gpu_mpi_main.cpp`,
  `faraday/ohm test_gpu_main.cpp` now covered by `tests/core/numerics/gpu`), the
  `boundary_condition/test_field_boundary_conditions_*.cpp`, `test_field_evolution.cpp`,
  `test_fields_flush_schedules.cpp`, `test_hybrid_level_initializer.cpp`, `test_geom.cpp` and
  `tests/diagnostic/test_diagnostics_{1d,2d}_{cmp,ref}.cpp`.
- Proper fix: register what still passes, delete what is superseded.

### Dead code
- `src/core/data/particles/mapper/` (bisection mappers, included nowhere, duplicates
  `CellFlattener`)
- `src/amr/solvers/threaded_solver_ppc.hpp` (609 lines) and `solver_ppc_threading_sorted.hpp`,
  `src/core/numerics/pusher/pusher.hpp`, included nowhere
- `.off` files: `magnetic_field_coarsener.hpp.off`, `interpolating.hpp.off`, `test_tagging.py.off`

### Personal paths in mkn.py
- `mkn.py`: `$HOME/git/phare/mkn`, `/home/deegan/intel/vtune/...`, `/opt/intel/oneapi/...`
- Proper fix: environment variables with defaults, or move it out of the repository.

### Pasted transcript below
- The SoAVX benchmark notes below the separator end in a pasted session transcript.
- Proper fix: keep the numbers and the conclusion, drop the transcript.

### 2D+ validations
- `doc/tiles.md` "what remains to be done": "2d+ validations (possible numeric issue)", not
  investigated yet.

##############################

SOAVX benchmark results


AVX doesn't pay off any more now that the PR's kernels are in. Here's SoAVX without and with AVX, in millions of particles per second:

┌────────────────┬────────┬──────┬─────────────────────────┐
│     SoAVX      │ no AVX │ AVX  │         Change          │
├────────────────┼────────┼──────┼─────────────────────────┤
│ 1D domain_only │ 26.4   │ 28.4 │ noise                   │
├────────────────┼────────┼──────┼─────────────────────────┤
│ 1D all         │ 28.9   │ 31.7 │ +9%, at the noise floor │
├────────────────┼────────┼──────┼─────────────────────────┤
│ 2D domain_only │ 21.0   │ 21.8 │ noise                   │
├────────────────┼────────┼──────┼─────────────────────────┤
│ 2D all         │ 26.8   │ 26.5 │ none                    │
├────────────────┼────────┼──────┼─────────────────────────┤
│ 3D domain_only │ 11.1   │ 10.5 │ noise                   │
├────────────────┼────────┼──────┼─────────────────────────┤
│ 3D all         │ 10.0   │ 7.8  │ −22%                    │
└────────────────┴────────┴──────┴─────────────────────────┘

What it tells us:
- The SIMD deposit is redundant now. In 2D, where it's the only AVX path running, it adds nothing. My SIMD deposit only vectorises the weight·coef·v products, and scatters with its own unfused per-corner indexing without __restrict. The PR's scalar fused deposit is now at least as fast.
- The 3D SIMD push is a clear regression (−22%), as before the PR, and now against a faster scalar baseline. m2p_avx builds large Array<Array<…>> temporaries and reads the fields with scattered loads, which the stride-based scalar m2p beats.
- The noise is large. Layouts AVX doesn't touch moved by up to 70% between the two runs (AoSTS 2D all went from 9.8 to 16.5). The bench warns that CPU scaling is enabled, and each configuration runs only 5 iterations, so treat anything under about 15% as noise.

My recommendation: drop the AVX paths from SoAVX for now:
- remove particleToMesh_avx from interpolating.hpp;
- make avx_ad_acc_ad always use the scalar loop, or have the backend call the scalar per-index loop directly.

SoAVX then runs the scalar fused path, which is what wins anyway. The pieces needed to compile with -w mkn.avx (the include and drift fixes, the loop-bound fix) can stay, since they're correct and harmless.

If you'd rather keep a SIMD deposit, the version worth trying is SIMD weights fed into the PR's fused deposit_ instead of my scatter. Given these numbers, I'd expect it to gain little.

To get firmer numbers, run with CPU scaling disabled (for example cpupower frequency-set -g performance) and more iterations (PHARE_PUSHES=20).

Which do you want: drop the AVX paths, or try SIMD weights fed into the fused deposit_?

✻ Worked for 13s · done 13:10

❯ ok so since https://github.com/PHAREHUB/PHARE/pull/1504 SoAVX is redundant?

