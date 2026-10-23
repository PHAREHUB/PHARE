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
