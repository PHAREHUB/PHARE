set -ex

# N_CORES should be set in docker agent as ncores available to agent
[ -z "$N_CORES" ] && echo "N_CORES not set: error" && exit 1

dnf install -y compiler-rt llvm elfutils # move to dockerfile
export ASAN_SYMBOLIZER_PATH=/usr/bin/llvm-symbolizer
export ASAN_OPTIONS=detect_leaks=0

GCC_ASAN_PRELOAD=$(gcc -print-file-name=libasan.so)
CLANG_ASAN_PRELOAD=$(clang -print-file-name=libclang_rt.asan.so)
cd build
ASAN_OPTIONS=detect_leaks=0
OMPI_MCA_memory=^patcher LD_PRELOAD=$CLANG_ASAN_PRELOAD ctest -j$N_CORES --output-on-failure
