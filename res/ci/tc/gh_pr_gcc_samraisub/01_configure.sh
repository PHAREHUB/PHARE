set -ex

mkdir -p build

cd build
cmake .. \
         -G Ninja -DCMAKE_CXX_FLAGS="-g3 -O3 -march=native -mtune=native -DPHARE_DIAG_DOUBLES=1" -DCMAKE_BUILD_TYPE=Debug
