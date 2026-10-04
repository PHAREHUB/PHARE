set -ex

# N_CORES should be set in docker agent as ncores available to agent
[ -z "$N_CORES" ] && echo "N_CORES not set: error" && exit 1

which mpirun

BUILD=$PWD

# teamcity-fedora_dep images have SAMRAI installed into the container
#  teamcity-fedora images do not, and SAMRAI is built as a subproject as needed
cmake . -DCMAKE_PREFIX_PATH=/root/.local-samrai -DexecuteNotebooks=ON \
        -Ddocumentation=ON \
        -DdevMode=ON -DtestMPI=ON \
        -DPHARE_EXEC_LEVEL_MAX=9 -Dphare_configurator=ON \
        -DCMAKE_CXX_FLAGS="-O3 -march=native"

ctest -j$(($N_CORES / 2)) --output-on-failure
