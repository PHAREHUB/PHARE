set -ex

mkdir -p build && cd build

# N_CORES should be set in docker agent as ncores available to agent
[ -z "$N_CORES" ] && echo "N_CORES not set: error" && exit 1
echo "building with N_CORES $N_CORES"
ninja -j10
