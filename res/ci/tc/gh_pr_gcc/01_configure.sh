set -ex

CMAKE_CXX_FLAGS="-O3 -DPHARE_DIAG_DOUBLES=1 "
cmake . -Ddocumentation=ON -Dphare_configurator=ON \
        -DdevMode=ON -G Ninja -DSAMRAI_ROOT=/usr/local \
        -DPHARE_EXEC_LEVEL_MAX=9 \
        -DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS}" \
        -DCMAKE_BUILD_TYPE=Debug
