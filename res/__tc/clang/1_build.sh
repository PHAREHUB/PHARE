set -ex

[ -z "$N_CORES" ] && echo "N_CORES not set: error" && exit 1

export CC=clang
export CXX=clang++
export ASAN_OPTIONS=detect_leaks=0

(
  mkdir build && cd build
  cmake .. -G Ninja -Dasan=ON  \
          -DdevMode=ON  \
          -DCMAKE_BUILD_TYPE=Debug \
          -Dphare_configurator=ON -DSAMRAI_ROOT=/usr/local \
          -DCMAKE_CXX_FLAGS="-g3 -O3 -march=native -mtune=native -DPHARE_DIAG_DOUBLES=1"

  # NO SUBPROJECT SAMRAI EXPECTED!
  [ -d ../subprojects/samrai ] && exit 1
  ninja -j$N_CORES -v
)
