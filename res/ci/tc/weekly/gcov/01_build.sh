set -ex

# N_CORES should be set in docker agent as ncores available to agent
[ -z "$N_CORES" ] && echo "N_CORES not set: error" && exit 1

git log -3
g++ -v

# no -DSAMRAI_ROOT: SAMRAI is built as a subproject
# no -G Ninja, cant use with coverage
# testMPI=OFF builds every test target, the MPI tests reuse them, see 02_test.sh
mkdir -p build && cd build
cmake .. -Dcoverage=ON -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_FLAGS="-O0" -DtestMPI=OFF
make -j$N_CORES
