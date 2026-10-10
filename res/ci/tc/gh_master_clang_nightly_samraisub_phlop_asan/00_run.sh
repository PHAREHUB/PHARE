set -ex

echo $N_CORES
[ -z "$N_CORES" ] && echo "N_CORES not set: error" && exit 1

python3 -m pip install psutil phlop vtk -U 
ulimit -n 66666

mkdir -p build

export OMPI_MCA_btl_vader_single_copy_mechanism=none

CMAKE_CXX_FLAGS="-O3 -g0 -DPHARE_DIAG_DOUBLES=1 -march=native -mtune=native -DNDEBUG "
CMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS} -DPHARE_FORCE_LOG_LINE=1 "
CMAKE_CONFIG=" -Dphare_configurator=ON "
CMAKE_BUILD_TYPE="Release"
CMAKE_CONFIG="${CMAKE_CONFIG} -DdevMode=ON -DSAMRAI_ROOT=/usr/local"

N_CORES=$(python3 -c "import psutil; print(psutil.cpu_count(logical=False))")
export CXX=clang++
export CC=clang
export ASAN_OPTIONS=detect_leaks=0
export CLANG_ASAN_PRELOAD=$(clang -print-file-name=libclang_rt.asan.so)
export OMPI_MCA_memory=^patcher 

(
  cd build && cmake .. ${CMAKE_CONFIG} -DtestMPI=ON -Dasan=ON \
       -DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS}" -G Ninja \
       -DCMAKE_BUILD_TYPE="${CMAKE_BUILD_TYPE}" -DPHARE_EXEC_LEVEL_MAX=100
  mold --run ninja -v -j${N_CORES}
)
python3 -m phlop.run.test_cases --cmake -i build --dump ctests.0.dill

(
  cd build 
  cmake .. ${CMAKE_CONFIG} -DtestMPI=OFF -Dasan=ON \
        -DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS}" -G Ninja \
        -DCMAKE_BUILD_TYPE="${CMAKE_BUILD_TYPE}" -DPHARE_EXEC_LEVEL_MAX=10
  mold --run ninja -v -j${N_CORES}
)

python3 -m phlop.run.test_cases --cmake -i build --dump ctests.1.dill
LD_PRELOAD=$CLANG_ASAN_PRELOAD python3 -m phlop.run.test_cases --load "ctests.*.dill" -Rc 32
