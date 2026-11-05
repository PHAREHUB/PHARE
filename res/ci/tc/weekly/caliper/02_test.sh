set -ex

# N_CORES should be set in docker agent as ncores available to agent
[ -z "$N_CORES" ] && echo "N_CORES not set: error" && exit 1

python3 -m pip install phlop -U

# per test caliper runtime report, kept in the .phlop logs
export CALI_CONFIG=runtime-report,calc.inclusive
python3 -m phlop.run.test_cases --cmake -i build --logging 2 -c $N_CORES
find .phlop
