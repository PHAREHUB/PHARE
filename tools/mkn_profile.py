#  running tests/core/numerics/ion_updater/test_multi_updater.cpp
#
#

import os
import sys
import shlex
import shutil
import subprocess
import numpy as np
from pathlib import Path


dim = 3
run_sim = 1  # switch to just plot and not run a sim
PATCHES = [1]
CELLS = [10]
PPC = [100]

permutables = [
    ("patches", PATCHES),
    ("cells", CELLS),
    ("ppc", PPC),
]

PHARE_ASYNC_TIMES = ".phare/async/multi_updater"

fn_strings = {
    0: "sync ghosts",
    1: "Boris::move_domain",
    2: "Boris::move_ghosts",
    3: "sync",
    4: "Group_barrier",
    5: "Domain_insert",
    6: "Group_barrier",
    7: "Deposit",
}

fn_strings = {
    0: "Boris::move_domain",
    1: "sync",
    2: "Domain_insert",
    3: "Deposit",
    # 4: "Reduce",
}


def run_permutation(patches, cells, ppc):
    times_dir = f"{PHARE_ASYNC_TIMES}/{patches}/{cells}/{ppc}"
    env = os.environ.copy()
    env_update = {
        "PHARE_ASYNC_TIMES": times_dir,
        "PHARE_PATCHES": f"{patches}",
        "PHARE_CELLS": f"{cells}",
        "PHARE_PPC": f"{ppc}",
    }
    env.update(env_update)
    print({k: v for k, v in env.items() if k.startswith("PHARE_")})

    run = "run" if not env.get("MKN_DBG") else "dbg"
    print("run", run, env.get("MKN_DBG"))

    if run_sim:
        cmd = "mkn run -p test_core"
        cmd = "python3 -O tests/functional/harris/harris_2d.py"
        subprocess.run(shlex.split(cmd), check=True, env=env)

    # shutil.copy(".phare_times.0.txt", f"{times_dir}/scope_times.txt")
    return times_dir


def parse_times(times_dir, patches, cells, ppc):
    times_per_type = {}

    for f in Path(times_dir).glob("*.txt"):
        type_id = f.stem
        if type_id == "scope_times":
            continue

        times = []
        with open(f) as file:
            for line in file.readlines():
                line = line.strip()
                times.append(line.split(" "))
        times_per_type[type_id] = times

    return {times_dir: (times_per_type, patches, cells, ppc)}


def plot_fn(times, fn):
    import matplotlib.pyplot as plt

    n_arrays = 2
    names = ["domain", "patchghost"]  # , "levelghost"]
    x_axis = [0.1, 0.2]  # , 0.2, 0.3]
    fig, ax = plt.subplots(figsize=(8.0, 8.0))

    for times_dir, times_per_type_tuple in times.items():
        times_per_type, patches, cells, ppc = times_per_type_tuple
        for type_id, bits_list in times_per_type.items():
            type_id_csv = type_id.split(",")
            type_str = type_id_csv[1]
            alloc_mode = type_id_csv[2][:3]
            type_str = f"{type_str}/{alloc_mode}/{patches}/{cells}/{ppc}"
            fn_0_times = [[], []]  # , []]

            pid = 0
            for bits in bits_list:
                if bits[1] == f"{fn}":
                    fn_times = fn_0_times[pid % n_arrays]
                    fn_times.append(int(bits[2]))
                    pid += 1

            print("\nfn_0_times", fn, type_id, fn_0_times)
            fn_0_times = [np.asarray(ts).mean() for ts in fn_0_times]
            ax.plot(x_axis, fn_0_times, "-", label=f"{type_str}")

    box = ax.get_position()
    ax.set_position(
        [box.x0, box.y0 + box.height * 0.2, box.width * 0.7, box.height * 0.8]
    )

    ax.legend(
        loc="center left",
        bbox_to_anchor=(1, 0.5),
        fancybox=True,
        shadow=True,
    )
    ax.set_title(fn_strings[fn])
    ax.set_xticks(x_axis)
    ax.set_xticklabels(names, rotation="vertical", fontsize=18)
    plt.ylabel("nanoseconds")
    fig.savefig(f"profile_plot.{fn}.png", dpi=100)


def plot(times):
    for i in range(len(fn_strings)):
        plot_fn(times, i)


def print_times(times):
    print("times")

    for _, info in times.items():
        nparts = (info[-2] ** dim) * info[-1]
        print(f"dim 3, cells {info[-2]}, ppc {info[-1]}, nparts {nparts}")
        for fn, stats in info[0].items():
            line_strs = [str(fn_strings[int(stat[1])]) for stat in stats]
            line_max = max([len(line_str) for line_str in line_strs])

            for i, line in enumerate(line_strs):
                padding = line_max - len(line)
                time = float(stats[i][2])
                norm_ppc = f"norm ppc: {time / nparts:>10,.2f}ns"
                print(f"{line}" + (" " * padding) + f"{time/1e6:>10,.2f}ms {norm_ppc}")


def main():
    import itertools

    times = {}
    permutations = itertools.product(*[e[1] for e in permutables])
    for permutation in permutations:
        times.update(parse_times(run_permutation(*permutation), *permutation))

    # plot(times)
    print_times(times)


if __name__ == "__main__":
    main()
