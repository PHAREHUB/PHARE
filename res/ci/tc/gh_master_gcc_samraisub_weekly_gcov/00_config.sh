set -ex

# N_CORES should be set in docker agent as ncores available to agent
[ -z "$N_CORES" ] && echo "N_CORES not set: error" && exit 1

git config --global --add safe.directory $PWD # weird git error resolution

mkdir -p build && cd build

# no -G Ninja, cant use with coverage
cmake .. -Dcoverage=ON \
         -DCMAKE_CXX_FLAGS="-O0" \
         -DCMAKE_BUILD_TYPE=Debug
