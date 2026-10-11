set -ex

# N_CORES should be set in docker agent as ncores available to agent
[ -z "$N_CORES" ] && echo "N_CORES not set: error" && exit 1

git log -3
g++ -v

# no -DSAMRAI_ROOT: SAMRAI is built as a subproject
# testMPI=OFF builds every test target, the MPI tests reuse them, see 02_test.sh
mkdir -p build && cd build
cmake .. -G Ninja -DdevMode=ON -Dphare_configurator=ON \
         -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_FLAGS="-O3 -DPHARE_DIAG_DOUBLES=1" \
         -DPHARE_EXEC_LEVEL_MAX=9 -DtestMPI=OFF
ninja -j$N_CORES
