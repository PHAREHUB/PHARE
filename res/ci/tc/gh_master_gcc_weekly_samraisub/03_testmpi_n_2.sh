set -ex

# N_CORES should be set in docker agent as ncores available to agent
[ -z "$N_CORES" ] && echo "N_CORES not set: error" && exit 1

mkdir -p build

CMAKE_CXX_FLAGS="-O3 -g0 -DPHARE_DIAG_DOUBLES=1 -march=native -mtune=native -DNDEBUG"
CMAKE_CONFIG=" -Dphare_configurator=ON -DPHARE_EXEC_LEVEL_MAX=100 " 
CMAKE_BUILD_TYPE="Release"
(
  cd build && cmake .. ${CMAKE_CONFIG} -DtestMPI=ON \
        -DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS}" -G Ninja \
        -DCMAKE_BUILD_TYPE="${CMAKE_BUILD_TYPE}" 
  ninja -v -j${N_CORES}
  ctest -j$((($N_CORES / 2) - 1)) --output-on-failure
)
