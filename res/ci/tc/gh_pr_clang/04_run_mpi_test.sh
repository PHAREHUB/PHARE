set -ex

# N_CORES should be set in docker agent as ncores available to agent
[ -z "$N_CORES" ] && echo "N_CORES not set: error" && exit 1

which mpirun

BUILD=$PWD

cd build
cmake .. -G Ninja -Dasan=ON -DtestMPI=ON \
        -DdevMode=ON -Dcppcheck=ON \
        -DCMAKE_CXX_FLAGS="-O3 -g3" -Dphare_configurator=ON

GCC_ASAN_PRELOAD=$(gcc -print-file-name=libasan.so)
CLANG_ASAN_PRELOAD=$(clang -print-file-name=libclang_rt.asan.so)
export ASAN_OPTIONS=detect_leaks=0
OMPI_MCA_memory=^patcher LD_PRELOAD=$CLANG_ASAN_PRELOAD ctest -j$(($N_CORES / 2)) --output-on-failure
