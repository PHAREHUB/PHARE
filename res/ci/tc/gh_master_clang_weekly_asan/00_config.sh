set -ex

# N_CORES should be set in docker agent as ncores available to agent
[ -z "$N_CORES" ] && echo "N_CORES not set: error" && exit 1

git config --global --add safe.directory $PWD # weird git error resolution
mkdir -p build && cd build

export CXX=clang++
export CC=clang
export ASAN_OPTIONS=detect_leaks=0

cmake .. -DtestMPI=ON -Dasan=ON -DPHARE_PYTEST_SIMULATORS=3 \
         -DCMAKE_CXX_FLAGS="-O3 -march=native -mtune=native -DPHARE_DIAG_DOUBLES=1 -DNDEBUG -g3" \
         -G Ninja -Dphare_configurator=ON \
         -DCMAKE_BUILD_TYPE=Release -DPHARE_EXEC_LEVEL_MAX=100
