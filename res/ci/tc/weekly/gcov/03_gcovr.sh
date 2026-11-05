set -ex

# coverage report from the data of 02_test.sh, into build/coverage
cd build
make gcovr
