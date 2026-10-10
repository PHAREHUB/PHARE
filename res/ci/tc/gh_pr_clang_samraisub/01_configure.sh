set -ex

mkdir -p build && cd build
echo $PWD

export CXX=clang++
export CC=clang

# C++20 module scanning needs clang-scan-deps, PHARE has no modules
cmake .. -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_SCAN_FOR_MODULES=0 \
         -DdevMode=ON \
         -G Ninja -DCMAKE_CXX_FLAGS="-O3 -DPHARE_DIAG_DOUBLES=1"
