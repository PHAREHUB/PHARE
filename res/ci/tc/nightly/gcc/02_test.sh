set -ex

# N_CORES should be set in docker agent as ncores available to agent
[ -z "$N_CORES" ] && echo "N_CORES not set: error" && exit 1

python3 -m pip install phlop psutil vtk -U
ulimit -n 66666
export OMPI_MCA_btl_vader_single_copy_mechanism=none

# both test lists dumped by 01_build.sh, tests needing the most cores first
python3 -m phlop.run.test_cases --load "ctests.*.dill" -Rc $N_CORES
