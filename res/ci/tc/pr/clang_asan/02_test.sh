set -ex

# N_CORES should be set in docker agent as ncores available to agent
[ -z "$N_CORES" ] && echo "N_CORES not set: error" && exit 1

dnf install -y compiler-rt llvm elfutils # move to dockerfile
python3 -m pip install phlop -U

export ASAN_SYMBOLIZER_PATH=/usr/bin/llvm-symbolizer
export ASAN_OPTIONS=detect_leaks=0
export OMPI_MCA_memory=^patcher

# serial tests, as built by 01_build.sh
python3 -m phlop.run.test_cases --cmake -i build --dump ctests.serial.dill

# MPI tests: only test registration changes with testMPI, the binaries are reused as is
(cd build && cmake .. -DtestMPI=ON)
python3 -m phlop.run.test_cases --cmake -i build --dump ctests.mpi.dill

# run both under the ASan runtime, tests needing the most cores first
CLANG_ASAN_PRELOAD=$(clang -print-file-name=libclang_rt.asan.so)
LD_PRELOAD=$CLANG_ASAN_PRELOAD python3 -m phlop.run.test_cases --load "ctests.*.dill" -Rc $N_CORES
