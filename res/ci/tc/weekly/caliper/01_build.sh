set -ex

# N_CORES should be set in docker agent as ncores available to agent
[ -z "$N_CORES" ] && echo "N_CORES not set: error" && exit 1

git log -3
g++ -v

# no -DSAMRAI_ROOT: SAMRAI is built as a subproject, -DwithCaliper downloads caliper too
# functional tests only, exec level 11 to 100
mkdir -p build && cd build
cmake .. -G Ninja -Dphare_configurator=ON -DwithCaliper=ON -DCMAKE_BUILD_TYPE=Release \
         -DCMAKE_CXX_FLAGS="-O3 -march=native -mtune=native -DPHARE_DIAG_DOUBLES=1 -DPHARE_LOG_LEVEL=1" \
         -DPHARE_EXEC_LEVEL_MIN=11 -DPHARE_EXEC_LEVEL_MAX=100
ninja -j$N_CORES
