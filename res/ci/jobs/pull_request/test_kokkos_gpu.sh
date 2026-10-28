#!/usr/bin/env bash
# GPU version of the "kokkos" job in .github/workflows/build.yml
#  the Kokkos backend is picked from the GPU found at runtime: CUDA (nvidia-smi) or HIP (rocminfo)
#  KOKKOS_ARCH=AMD_GFX1100 (for example) overrides the detected architecture

(
  set -ex

  CWD="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd)"
  BUILD_DIR=${BUILD_DIR:-"${CWD}/build"}

  if command -v nvidia-smi >/dev/null && nvidia-smi -L >/dev/null 2>&1; then
    BACKEND="-DKokkos_ENABLE_CUDA=ON"
    export CC=${CC:-clang} CXX=${CXX:-clang++}
  elif command -v rocminfo >/dev/null && rocminfo | grep -q "gfx"; then
    GFX=$(rocminfo | grep -o -m1 "gfx[0-9a-f]*")
    KOKKOS_ARCH=${KOKKOS_ARCH:-"AMD_${GFX^^}"}
    BACKEND="-DKokkos_ENABLE_HIP=ON"
    export CC=${CC:-hipcc} CXX=${CXX:-hipcc}
    export ROCM_PATH=${ROCM_PATH:-$(hipconfig --rocmpath)}
  else
    echo "no CUDA or HIP GPU found" >&2
    exit 1
  fi
  [ -n "${KOKKOS_ARCH}" ] && BACKEND="${BACKEND} -DKokkos_ARCH_${KOKKOS_ARCH}=ON"

  [ -f /.venv/bin/activate ] && . /.venv/bin/activate

  mkdir -p "${BUILD_DIR}" && cd "${BUILD_DIR}"
  cmake "${CWD}" -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON --fresh \
        -DCMAKE_BUILD_TYPE=Debug -DdevMode=ON -DwithKokkos=ON -DwithGpu=ON ${BACKEND} \
        -DCMAKE_CXX_FLAGS="-O2 -DPHARE_DIAG_DOUBLES=1" ${SAMRAI_ROOT:+-DSAMRAI_ROOT=${SAMRAI_ROOT}}
  cmake --build . --target test-field-solvers-gpu -j ${THREADS:-$(nproc)}
  ctest -R test-field-solvers-gpu --output-on-failure
)
