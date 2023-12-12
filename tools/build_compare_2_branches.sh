
set -ex
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $SCRIPT_DIR && cd .. && CWD=$PWD # move to project root

python3 tools/build_compare_2_branches.py $@
