set -ex

# N_CORES should be set in docker agent as ncores available to agent
[ -z "$N_CORES" ] && echo "N_CORES not set: error" && exit 1

mkdir -p build && cd build

cmake .. \
         -DCMAKE_CXX_FLAGS="-O3 -march=native -mtune=native" \
         -DtestMPI=ON \
         -G Ninja 
ctest -j$((($N_CORES / 2) - 1)) --output-on-failure
