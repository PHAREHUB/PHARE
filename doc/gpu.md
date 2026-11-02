# GPU notes

How the GPU builds are set up and the rules for writing device code in PHARE. Open work and known
shortcuts are in [todo_gpu.md](todo_gpu.md). Paths below are examples from one machine.

## Builds

| build | how | what it is |
|---|---|---|
| CMake, real GPU | `tools/cmake.sh`, `CC=CXX=hipcc`, `-DwithMkn=ON -DwithGpu=ON` | `PHARE_HAVE_GPU_HARDWARE=1`, mkn.gpu kernels on the device |
| CMake, pseudo GPU | `-DwithMkn=ON` without `withGpu` | `MKN_GPU_CPU=1`: `GPU_UNIFIED` code runs on mkn.gpu's CPU backend |
| mkn, real GPU | `./mkn_rocm.sh -Og 0` (`res/mkn/hip.yaml`) | builds and runs `FILE` from the script |
| mkn, pseudo GPU | `./mkn.sh` (`res/mkn/mpi.yaml`, `MKN_GPU_CPU=1`, `PHARE_HAVE_THRUST=0`) | same tests on the CPU backend |
| Kokkos GPU | `res/ci/jobs/pull_request/test_kokkos_gpu.sh` | `-DwithKokkos=ON -DwithGpu=ON`, backend from the GPU found (CUDA or HIP) |

- `withGpu` needs `withMkn` or `withKokkos` (`res/cmake/def.cmake`). Kokkos alone can be CPU only.
- The Kokkos job builds Kokkos from source. Set `SAMRAI_ROOT` to an installed SAMRAI (else it is
  built too, with the GPU compiler) and `BUILD_DIR` to keep it out of `build/`. KokkosTools' roctx
  connector needs `ROCM_PATH` pointing at a ROCm with roctx (e.g. `/opt/rocm-6.2.1`).
- mkn builds both write `bin/core/libphare_core.a` and `mkn clean` wipes `bin/core/`. A library left
  by the other build fails to link with `-fPIC` errors: delete it to rebuild.
- Fedora ROCm puts `libamdhip64` in `/usr/lib64`, clang's HIP auto link then drops `-lamdhip64`:
  it is given explicitly in `res/mkn/hip*.yaml`, and the local (git ignored) `tools/cmake.sh` needs
  `-DCMAKE_EXE_LINKER_FLAGS=-lamdhip64 -DCMAKE_SHARED_LINKER_FLAGS=-lamdhip64`.
- Thrust (rocThrust) is auto detected (`core/def/detail/thrust.hpp`), some GPU particle services
  need it, so the pseudo GPU build (`PHARE_HAVE_THRUST=0`) can not run everything.

### mkn.gpu CPU backend (pseudo GPU)

Since mkn.gpu `b06035b` (#37) the CPU backend walks the full 3D grid and block, sets thread local
`threadIdx`/`blockIdx`/`blockDim`/`gridDim` (one host thread per stream runs kernels concurrently)
and has serial `atomicAdd`/`atomicSub`/`atomicCAS`. Warp size is 1 there, so PHARE launches one
thread per block and the no-op `__syncthreads` is correct. There is no dynamic shared memory:
kernels using `extern __shared__` stay behind `PHARE_HAVE_MKN_GPU_HW`.

## Layouts

- AoSPCTS is the `GPU_UNIFIED` layout under test.
- AoSCMTS is CPU only (its per tile cell map lives on the host).
- AoSTS is a proof of concept with the fused push + deposit kernel, not exposed to Python or tests.
- `IonUpdaterFor_t` picks `mkn_xyz::ParallelIonUpdater` for `GPU_UNIFIED` (not AoSMapped) when
  mkn.gpu is available, the host updaters otherwise.

## Writing device code

- Guard mkn stream kernels on `PHARE_HAVE_MKN_GPU` (hardware or pseudo GPU), not
  `PHARE_HAVE_MKN_GPU_HW`, unless they need shared memory.
- `PHARE_DEVICE_PASS` (`core/def/phare_config.hpp`) is 1 only in the device compilation pass:
  branch on it for device builtins (`atomicAdd`, `warpSize`) and for host only code in
  `_PHARE_ALL_FN_` functions. `__HIPCC__`/`__CUDACC__` are set in both passes.
- `throw_exception` asserts in the device pass and throws on the host. Device code reports errors
  through a record checked on the host afterwards, e.g. `boris::MoveError` and
  `MultiBoris::check_move_error`.
- Kernel lambdas are `[=] _PHARE_ALL_FN_(...)`, the attribute before the parameter list. Capture
  copies of views, never `[&]` or `core::Field` (it holds a `std::string`): decay with CTAD,
  `basic::Field f = field;`, and for vector fields build a `basic::TensorField` of those.
- Copy `*this` explicitly (`self = *this`) when a kernel needs the solver/functor.
- In a function with a trailing `requires`, the attribute goes after it:
  `auto f(Args&&...) requires(...) _PHARE_ALL_FN_`.
- `constexpr` functions are implicitly host + device with clang; deduction guides are not, e.g. use
  `std::array<T, N>{...}`, not `std::array{...}`.
- Gate host only debug code on the instantiation's own trait (`if constexpr (alloc_mode != ...)`)
  when the type exposes it, else remove it rather than adding trait plumbing.
- Kokkos `parallel_for` is asynchronous: `Kokkos::fence()` before the host reads or writes the data.
- Results may differ from the CPU by a few ulps (FMA contraction on the device): compare with a
  tolerance, relative to the field's magnitude when it is large.

## Tests

- Backend agnostic: compare the CPU result with `GPU_UNIFIED`, register `if(withGpu)`, no Kokkos,
  mkn or RAJA specific tests. See `tests/core/numerics/gpu/test_field_solvers.cpp`.
- Small targeted tests for one mechanism, e.g. `tests/core/numerics/pusher/test_boris_move_error.cpp`,
  rather than more cases in `test_multi_updater.cpp`.
- `test_multi_updater.cpp` has AoSPCTS `GPU_UNIFIED` cases under `PHARE_WITH_MKN_GPU`, so they run
  on the pseudo GPU too.

## Shortcuts

Every GPU shortcut (host fallback, fixed default, excluded layout) gets an entry in
[todo_gpu.md](todo_gpu.md): where, what, and the proper fix, in the change that introduces it.
