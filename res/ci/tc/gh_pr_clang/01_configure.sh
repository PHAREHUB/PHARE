set -ex
pwd
ulimit -a

export CXX=clang++
export CC=clang
export ASAN_OPTIONS=detect_leaks=0

mkdir build && cd build
# C++20 module scanning needs clang-scan-deps, PHARE has no modules
cmake .. -G Ninja -Dasan=ON -DCMAKE_CXX_SCAN_FOR_MODULES=0 \
        -DdevMode=ON -Dcppcheck=ON -DSAMRAI_ROOT=/usr/local \
        -DCMAKE_CXX_FLAGS="-O3 -g3 -DPHARE_DIAG_DOUBLES=1" -DCMAKE_BUILD_TYPE=Debug -Dphare_configurator=ON
cd ..
unset LD_PRELOAD
# should be system version                
[ -d subprojects/samrai ] && exit 1 

exit 0
