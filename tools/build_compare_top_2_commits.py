import shutil
from pathlib import Path

from pyphare.simulator.simulators import Simulators
from tools.python3 import pushd, cmake, git
from tests.functional.alfven_wave import alfven_wave1d

# we want it seeded
alfven_wave1d.MODEL_INIT={"seed": 1337}
alfven_wave1d.TIME_STEP_NBR = 10

# register exit handler
git.git_branch_reset_at_exit()

top_2_hashes  = git.hashes(2)

build_dir = Path("build")

for hsh in top_2_hashes:
    b = (build_dir / hsh)
    b.mkdir(parents=True, exist_ok=True)

    git.checkout(hsh)
    with pushd(b):
        cmake.config("../..")

for hsh in top_2_hashes:
    b = (build_dir / hsh)

    git.checkout(hsh)
    with pushd(b):
        cmake.build()

simulators = Simulators()
for hsh in top_2_hashes:
    simulators.register(alfven_wave1d.config(), build_dir=str(build_dir / hsh))

simulators.run(compare=True)
