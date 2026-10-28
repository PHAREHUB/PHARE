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

## Tiling

### 3D GPU tile size is a fixed default
- `src/core/data/tiles/tiling_options.hpp` (`TilingOptions::defaults_for`)
- `GPU_UNIFIED` in 3D defaults to `max_tile_size = 6` so a tile's ghost box of field nodes fits
  the GPU kernels' shared memory (`ds * 4 * 8` + warp buffer, asserted < 65000). Sized for
  interp order 1 only.
- Proper fix: derive the max tile size from the ghost width (interp order) and the device's
  shared memory limit.

## Boundary conditions

### Tiled physical boundary conditions gather/scatter a whole patch
- `src/amr/data/field/refine/field_refine_patch_strategy.hpp` (`setPhysicalBoundaryConditions`)
- Tiled fields are reduced into a temporary untiled grid, the boundary condition is applied there,
  and the result is scattered back into every tile. Allocates per call and is host only.
- Proper fix: apply boundary conditions per tile on the tile's ghost region.

## Field solvers

### Faraday and Ohm have no device kernels
- `src/core/numerics/faraday/faraday.hpp`, `src/core/numerics/ohm/ohm.hpp`
- Their kernels are `[&]` capturing host lambdas over master's `curl(layout_, E)` /
  `grad(layout_, Pe)` (#1468), so `GPU_UNIFIED` field solves run on the host. Ampere is done:
  `curl` is `_PHARE_ALL_FN_` and its lambdas copy `basic::Field` views of J.
- Proper fix: same as Ampere, copy `basic::Field` views into `[=] _PHARE_ALL_FN_` lambdas, and make
  `grad` and Ohm's `ideal_`/`pressure_`/... device callable.



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

