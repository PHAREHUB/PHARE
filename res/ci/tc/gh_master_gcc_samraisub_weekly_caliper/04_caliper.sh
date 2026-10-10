set -ex

# N_CORES should be set in docker agent as ncores available to agent
[ -z "$N_CORES" ] && echo "N_CORES not set: error" && exit 1

rm -rf build && mkdir build && cd build

cmake .. \
         -DCMAKE_CXX_FLAGS="-O3 -march=native -mtune=native" \
         -G Ninja -DwithCaliper=1 \
         -DCMAKE_BUILD_TYPE=Release

ninja
ctest -j${N_CORES}
export CALI_CONFIG=runtime-report
ctest -R py3_initialization -V

cd ..
export PYTHONPATH=$PWD/build:$PWD/pyphare
mpirun -n 2 python3 -u tests/simulator/test_initialization.py
mpirun -n 3 python3 -u tests/simulator/test_initialization.py
mpirun -n 4 python3 -u tests/simulator/test_initialization.py
mpirun -n 5 python3 -u tests/simulator/test_initialization.py
