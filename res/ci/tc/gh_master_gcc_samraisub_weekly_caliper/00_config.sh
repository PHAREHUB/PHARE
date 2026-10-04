set -ex

git config --global --add safe.directory $PWD # weird git error resolution

git log -3

ulimit -a

# N_CORES should be set in docker agent as ncores available to agent
[ -z "$N_CORES" ] && echo "N_CORES not set: error" && exit 1

mkdir -p build 

(cd build && cmake .. -Dphare_configurator=ON \
         -DCMAKE_CXX_FLAGS="-O3 -march=native -mtune=native -DPHARE_DIAG_DOUBLES=1 -DPHARE_LOG_LEVEL=1" \
         -G Ninja -DwithCaliper=ON \
         -DCMAKE_BUILD_TYPE=Release -DPHARE_EXEC_LEVEL_MIN=11 -DPHARE_EXEC_LEVEL_MAX=100)
