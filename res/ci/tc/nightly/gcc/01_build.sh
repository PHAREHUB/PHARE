set -ex

# N_CORES should be set in docker agent as ncores available to agent
[ -z "$N_CORES" ] && echo "N_CORES not set: error" && exit 1

git log -3
g++ -v

python3 -m pip install phlop -U

# no -DSAMRAI_ROOT: SAMRAI is built as a subproject
mkdir -p build && cd build

# MPI tests up to exec level 100, including the functional tests
cmake .. -G Ninja -DdevMode=ON -Dphare_configurator=ON -DCMAKE_BUILD_TYPE=Release \
         -DCMAKE_CXX_FLAGS="-O3 -g0 -march=native -mtune=native -DNDEBUG -DPHARE_DIAG_DOUBLES=1 -DPHARE_FORCE_LOG_LINE=1" \
         -DtestMPI=ON -DPHARE_EXEC_LEVEL_MAX=100
ninja -j$N_CORES
python3 -m phlop.run.test_cases --cmake -i . --dump ../ctests.mpi.dill

# serial tests up to exec level 10, same flags (cached)
cmake .. -DtestMPI=OFF -DPHARE_EXEC_LEVEL_MAX=10
ninja -j$N_CORES
python3 -m phlop.run.test_cases --cmake -i . --dump ../ctests.serial.dill
