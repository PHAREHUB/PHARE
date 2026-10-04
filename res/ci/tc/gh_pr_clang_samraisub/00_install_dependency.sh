set -ex
# this is for TC errors
git config --global --add safe.directory $PWD
clang++ -v
