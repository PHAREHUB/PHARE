
# our tests can require files from github/etc

import os, subprocess

dir_path = os.path.dirname(os.path.realpath(__file__))
overlaps_data = os.path.join(dir_path, "overlaps_test_data")

def run_shell(cmd):
    return subprocess.run(cmd, shell=True, check=True)

if not os.path.isdir(overlaps_data):
    run_shell("git clone https://github.com/PHARCHIVE/hdf5_patch_overlap --depth 1 overlaps_test_data")
    run_shell("find . -not -path '*/\\.*' -type f -name '*.xz' | xargs unxz")
