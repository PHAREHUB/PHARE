set -ex

[ -z "$N_CORES" ] && echo "N_CORES not set: error" && exit 1

(
  mkdir build && cd build
  cmake .. -G Ninja  \
          -DdevMode=ON  \
          -DCMAKE_BUILD_TYPE=Debug \
          -Dphare_configurator=ON \
          -DCMAKE_CXX_FLAGS="-g3 -O3 -march=native -mtune=native -DPHARE_DIAG_DOUBLES=1"

  ninja -j$N_CORES -v
)
