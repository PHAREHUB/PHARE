set -ex

cd samrai/build
# configure script expect to find unordered_map in tr1
mkdir /usr/include/boost/tr1
ln -s /usr/include/boost/unordered_map.hpp /usr/include/boost/tr1

sh ../SAMRAI/configure --with-boost=/usr --with-CXX=mpicxx --with-CC=mpicc --with-F77=mpif77 --with-hdf5=/usr
make dox
