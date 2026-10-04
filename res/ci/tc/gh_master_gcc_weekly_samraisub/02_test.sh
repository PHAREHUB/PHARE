set -ex

mkdir -p build && cd build
ctest -j$N_CORES --output-on-failure
