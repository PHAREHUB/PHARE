set -ex

# N_CORES should be set in docker agent as ncores available to agent
[ -z "$N_CORES" ] && echo "N_CORES not set: error" && exit 1

git log -3
clang++ -v

export CXX=clang++
export CC=clang

# no -DSAMRAI_ROOT: SAMRAI is built as a subproject
# testMPI=OFF builds every test target, the MPI tests reuse them, see 02_test.sh
# C++20 module scanning needs clang-scan-deps, PHARE has no modules
mkdir -p build && cd build
cmake .. -G Ninja -Dasan=ON -DCMAKE_CXX_SCAN_FOR_MODULES=0 -DdevMode=ON -Dphare_configurator=ON \
         -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_FLAGS="-O3 -g3 -DPHARE_DIAG_DOUBLES=1" \
         -DPHARE_EXEC_LEVEL_MAX=9 -DtestMPI=OFF
ninja -j$N_CORES
