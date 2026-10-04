set -ex

# N_CORES should be set in docker agent as ncores available to agent
[ -z "$N_CORES" ] && echo "N_CORES not set: error" && exit 1

cd build

ctest -j$N_CORES --output-on-failure

BUILD=$PWD

cmake .. -Dcoverage=ON -DtestMPI=ON \
         -DCMAKE_CXX_FLAGS="-O0" \
         -DCMAKE_BUILD_TYPE=Debug

ctest -j$((($N_CORES / 2) - 1)) --output-on-failure

find tests/simulator -type f
