set -ex
pwd
ulimit -a

[ -z "$N_CORES" ] && echo "N_CORES not set: error" && exit 1

export TMPDIR=/tmp
export MODULEPATH=/etc/scl/modulefiles:/etc/scl/modulefiles:/usr/share/Modules/modulefiles:/etc/modulefiles:/usr/share/modulefiles
export GTEST_OUTPUT=xml:gtest_out.xml

# prevent blas_thread_server from spawning a bajillion threads
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1


export ASAN_SYMBOLIZER_PATH=/usr/bin/llvm-symbolizer
export ASAN_OPTIONS=detect_leaks=0
GCC_ASAN_PRELOAD=$(gcc -print-file-name=libasan.so)
CLANG_ASAN_PRELOAD=$(clang -print-file-name=libclang_rt.asan.so)

(
  cd build
  cmake .. -G Ninja -Dasan=ON -DtestMPI=ON \
          -DdevMode=ON -Dcppcheck=ON \
          -DCMAKE_CXX_FLAGS="-O3 -g3" -Dphare_configurator=ON
  ninja -j$N_CORES -v
  LD_PRELOAD=$CLANG_ASAN_PRELOAD ctest -j$N_CORES --output-on-failure
)
