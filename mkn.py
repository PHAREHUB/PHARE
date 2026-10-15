#!/usr/bin/env python3
"""
export PYTHONPATH="$HOME/git/phare/mkn:$HOME/git/phare/mkn/build:$HOME/git/phare/mkn/pyphare:$HOME/git/phlop"
./mkn.sh  -Og 0  clean build -w mkn.gpu  -a -fno-omit-frame-pointer
./tools/mkn_profile.sh
python3 tools/python3/phloping.py print_scope_timings -f .phare_times.0.txt
"""


import os
import shutil
import subprocess

compile = "./mkn.sh -Og 0 clean build -w mkn.gpu -a -fno-omit-frame-pointer run"
procs = [
    "./tools/mkn_profile.sh",
    "python3 tools/python3/phloping.py print_scope_timings -f .phare_times.0.txt",
]

record = "perf record -e cycles,instructions,cache-misses,L1-dcache-load-misses -Tga -F 10000"
_update = {"MKN_DBG": record}


def run(proc, env=os.environ.copy()):
    subprocess.run(proc.split(" "), check=True, env=env)


def do(env_update):
    env = os.environ.copy()
    env.update(_update)
    env.update(env_update)
    for proc in procs:
        run(proc, env)


def vtune(name):
    dir = f"/home/deegan/intel/vtune/projects/p/{name}"
    shutil.rmtree(dir)
    run(
        f"/opt/intel/oneapi/vtune/latest/bin64/vtune -import {name}.perf -result-dir {dir}"
    )


def ref():
    do({"MKN_DBG": record + " -o ref.perf", "PHARE_REF_ONLY": "1"})
    # vtune("ref")


def cmp():
    do({"MKN_DBG": record + " -o cmp.perf", "PHARE_CMP_ONLY": "1"})
    # vtune("cmp")


# /opt/intel/oneapi/vtune/latest/bin64/vtune -import  cmp.perf -result-dir   /home/deegan/intel/vtune/projects/p/cmp
def main():
    run(compile)
    ref()
    cmp()


if __name__ == "__main__":
    main()
