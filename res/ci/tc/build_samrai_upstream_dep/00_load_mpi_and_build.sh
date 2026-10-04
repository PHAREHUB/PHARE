set -ex

# This build is purely for SAMRAI documentation creation

doxygen -v

GIT_URL="https://github.com/llnl/samrai"
VERSION="develop"
git clone --depth 1 $GIT_URL -b $VERSION samrai --recursive 
cd samrai

mkdir -p build && cd build

      
CMAKE_CXX_FLAGS="-g3 -O0 -fno-omit-frame-pointer"
CMAKE_BUILD_TYPE="Debug"
CMAKE_CONFIG+=" -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}"
      
cmake .. \
      -DENABLE_OPENMP=OFF -DENABLE_SAMRAI_TESTS=OFF \
      -DENABLE_DOCS=ON \
      -DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=true \
      -DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS}" \
      ${CMAKE_CONFIG}

[ -z "$N_CORES" ] && echo "N_CORES not set: error" && exit 1
echo "building with N_CORES $N_CORES"
make VERBOSE=1 -j$((N_CORES + 2)) && make install
make doxygen_docs
