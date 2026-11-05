set -ex

# functional tests (exec level 11+) with PHARE_DRY_RUN=1: set up, don't advance
./res/ci/jobs/pull_request/dry_run_functional_tests.sh
