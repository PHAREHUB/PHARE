set -ex

mkdir -p build && cd build
echo $PWD
ls -l

ninja -j10
