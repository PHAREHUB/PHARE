set -ex

[ -z "$N_CORES" ] && echo "N_CORES not set: error" && exit 1

git config --global --add safe.directory $PWD # weird git error resolution

export CALI_CONFIG=runtime-report,calc.inclusive

python3 -m pip install phlop -U
python3 -m phlop.run.test_cases --cmake -i build --logging 2 -c 32
find .phlop
