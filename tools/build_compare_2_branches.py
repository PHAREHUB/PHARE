import sys
import shutil
from pathlib import Path

from pyphare.simulator.simulators import Simulators
from tools.python3 import pushd, cmake, git
from tests.functional.alfven_wave import alfven_wave1d

# we want it seeded
alfven_wave1d.MODEL_INIT={"seed": 1337}
alfven_wave1d.TIME_STEP_NBR = 10

if len(sys.argv) != 3:
    print("Incorrect input arguments, expects two branch names", sys.argv)
    exit(0)

# register exit handler
git.git_branch_reset_at_exit()

branches  = sys.argv[1:]

# check branches exist
for branch in branches:
    git.checkout(branch)

build_dir = Path("build")

for branch in branches:
    b = (build_dir / branch)
    b.mkdir(parents=True, exist_ok=True)

    git.checkout(branch)
    with pushd(b):
        cmake.config("../..")

for branch in branches:
    git.checkout(branch)
    with pushd(build_dir / branch):
        cmake.build()

simulators = Simulators()
for branch in branches:
    # alfven_wave1d.config() will be identical here even if different on branches
    #  as it is already parsed before we change branch
    simulators.register(alfven_wave1d.config(), build_dir=str(build_dir / branch))
simulators.run(compare=True)
