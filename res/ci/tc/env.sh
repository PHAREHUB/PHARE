# sourced by every TeamCity build step before the step script, see <build>.kt

export TMPDIR=/tmp
export GTEST_OUTPUT=xml:gtest_out.xml

# prevent blas_thread_server from spawning a bajillion threads
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# docker agents run as root, sudo mpirun is not advised otherwise
export OMPI_ALLOW_RUN_AS_ROOT=1
export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

export MODULEPATH=/etc/scl/modulefiles:/usr/share/Modules/modulefiles:/etc/modulefiles:/usr/share/modulefiles
if command -v modulecmd > /dev/null; then
  eval `modulecmd bash load mpi/openmpi-x86_64`
fi
