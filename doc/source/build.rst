
===========
Build PHARE
===========



Build for production
--------------------

.. code-block::bash

        cd path/to/dir/containing/PHARE
        mkdir build
        cd build
        cmake -DCMAKE_CXX_FLAGS="-O3 -march=native -mtune=native" -DCMAKE_BUILD_TYPE=Release   ../PHARE
        make -j


Build for debugging
-------------------

.. code-block::bash

        cd path/to/dir/containing/PHARE
        mkdir build
        cd build
        cmake -DCMAKE_CXX_FLAGS="-g3 -O0 -march=native -mtune=native" -DCMAKE_BUILD_TYPE=Debug  -DCMAKE_CXX_FLAGS="-DPHARE_DIAG_DOUBLES=1"  ../PHARE
        make -j
