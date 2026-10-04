set -ex

ROOT=$PWD

cd build

# N_CORES should be set in docker agent as ncores available to agent
[ -z "$N_CORES" ] && echo "N_CORES not set: error" && exit 1

ctest -j$N_CORES --output-on-failure

cd ..

DRY_RUN="res/ci/jobs/pull_request/dry_run_functional_tests.sh"
if [ -f "${DRY_RUN}" ]; then
  "./${DRY_RUN}"
fi
