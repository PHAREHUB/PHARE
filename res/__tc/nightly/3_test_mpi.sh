set -ex

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


(
  cd build
  cmake .. -G Ninja -DtestMPI=ON \
          -DdevMode=ON -Dcppcheck=ON \
          -DCMAKE_CXX_FLAGS="-O3 -g3" -Dphare_configurator=ON
  ninja -j$N_CORES -v
  ctest -j$N_CORES --output-on-failure
)
