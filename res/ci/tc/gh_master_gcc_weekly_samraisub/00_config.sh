set -ex

# N_CORES should be set in docker agent as ncores available to agent
[ -z "$N_CORES" ] && echo "N_CORES not set: error" && exit 1

git config --global --add safe.directory $PWD # weird git error resolution
mkdir -p build && cd build

cmake .. -DPHARE_PYTEST_SIMULATORS=3 \
         -DCMAKE_CXX_FLAGS="-O3 -march=native -mtune=native -DPHARE_DIAG_DOUBLES=1 -DNDEBUG -g0" \
         -G Ninja -Dphare_configurator=ON \
         -DCMAKE_BUILD_TYPE=Release -DPHARE_EXEC_LEVEL_MAX=100
