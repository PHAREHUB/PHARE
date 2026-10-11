set -ex

# N_CORES should be set in docker agent as ncores available to agent
[ -z "$N_CORES" ] && echo "N_CORES not set: error" && exit 1

python3 -m pip install phlop -U

# serial tests, as built by 01_build.sh
python3 -m phlop.run.test_cases --cmake -i build --dump ctests.serial.dill

# MPI tests: only test registration changes with testMPI, the binaries are reused as is
(cd build && cmake .. -DtestMPI=ON)
python3 -m phlop.run.test_cases --cmake -i build --dump ctests.mpi.dill

# run both, tests needing the most cores first
python3 -m phlop.run.test_cases --load "ctests.*.dill" -Rc $N_CORES
