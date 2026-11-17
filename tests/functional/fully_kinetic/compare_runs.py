#!/usr/bin/env python3
"""
Fully kinetic static refined box (the configuration of tests/simulator/test_fully_kinetic_amr.py,
T2), run twice and compared on both levels (P4 tests T3c, T3d):

    --tiling        one patch per level (largest_patch_size 32) vs tiled (8): fine patch corners
                    on the level boundary, so the level ghost band and the same-level border sums
                    meet there
    --ranks N       tiled, serial vs N ranks (mpirun)
    --compare A B   two existing diagnostic directories

The arguments after "--" are T2's: dimension, interpolation order, --filter-passes, --cleaning.
E and B are compared bit for bit on the physical nodes of every patch of both levels at every
dumped time (single precision dumps), the domain particles of both levels sorted, counts exactly
and values within --particle-tolerance (the border sums add in an order that depends on the
tiling: round-off). Exit code 0 when everything matches.

Measured at d78ea99a8 (Debug), 2D, s = 0.5, 20 coarse steps:
    order 2, no filter: tiling and 3 ranks identical, particles within 8e-15
    order 1, or with the filter: differ, the level ghost band holds the particles of one cell
    beyond the patch only at order 2 (P4 plan, band widening)

    pbuild.sh run python3 -u tests/functional/fully_kinetic/compare_runs.py --tiling -- 2 2 \
        --filter-passes 0 --cleaning 0.5
"""

import argparse
import itertools
import os
import subprocess
import sys

import h5py
import numpy as np

ONE_PATCH, TILED = 32, 8


def run(patch_size, diag_dir, t2_args):
    """one T2 simulation with the given largest patch size, dumps in diag_dir"""
    sys.argv = [sys.argv[0]] + t2_args
    from tests.simulator import test_fully_kinetic_amr as amr
    from pyphare.simulator.simulator import Simulator
    import pyphare.pharein as ph

    amr.fk.CASES[amr.ndim]["tiled"] = patch_size
    amr.fk.diag_dir = diag_dir
    ph.global_vars.sim = None
    Simulator(amr.config()).run().reset()


def launch(patch_size, diag_dir, t2_args, ranks=1):
    """the same, in a fresh process (one simulation per process), on ranks MPI ranks"""
    command = [sys.executable, "-u", os.path.abspath(__file__), "--run", str(patch_size),
               diag_dir, "--", *t2_args]
    if ranks > 1:
        command = ["mpirun", "-n", str(ranks)] + command
    subprocess.run(command, check=True)


def dump_times(path):
    with h5py.File(f"{path}/EM_E.h5", "r") as f:
        return sorted(float(t) for t in f["t"])


def physical_nodes(path, quantity, time, level):
    """{(component, 2 x / dl, ...): value} on the physical nodes of every patch of the level"""
    from pyphare.pharesee.hierarchy import hierarchy_from

    hier = hierarchy_from(h5_filename=f"{path}/EM_{quantity}.h5", times=f"{time:.10f}")
    values = {}
    if level not in hier.levels(time):
        return values
    for patch in hier.level(level, time).patches:
        dl = patch.layout.dl
        lower = [lo * d for lo, d in zip(patch.box.lower, dl)]
        upper = [(up + 1) * d for up, d in zip(patch.box.upper, dl)]
        for component in "xyz":
            pdata = patch.patch_datas[f"{quantity}{component}"]
            axes = [getattr(pdata, axis) for axis in "xyz"[: len(dl)]]
            data = pdata.dataset[:]
            for index in itertools.product(*[range(len(axis)) for axis in axes]):
                position = [axis[i] for axis, i in zip(axes, index)]
                if all(lo - 1e-9 <= p <= up + 1e-9 for p, lo, up in zip(position, lower, upper)):
                    key = (component,) + tuple(round(2 * p / d) for p, d in zip(position, dl))
                    values[key] = data[index]
    return values


def sorted_particles(path, population, time_key, level):
    """iCell, delta, v, weight of the domain particles of the level, rows sorted"""
    rows = []
    with h5py.File(f"{path}/ions_pop_{population}_domain.h5", "r") as f:
        group = f["t"][time_key].get(f"pl{level}")
        for patch in group.values() if group is not None else []:
            n = patch["weight"][()].reshape(-1).size
            rows.append(
                np.hstack(
                    [patch[k][()].reshape(n, -1).astype(np.float64)
                     for k in ["iCell", "delta", "v", "weight"]]
                )
            )
    if not rows:
        return np.zeros((0, 0))
    array = np.vstack(rows)
    return array[np.lexsort(array.T[::-1])]


def compare(a, b, particle_tolerance):
    """True when E, B are bitwise equal and the particles match, prints a summary"""
    ok = True
    times = dump_times(a)
    if times != dump_times(b):
        print("dump times differ")
        return False

    values = differ = 0
    first = None
    for time in times:
        for level in [0, 1]:
            for quantity in ["E", "B"]:
                va, vb = physical_nodes(a, quantity, time, level), physical_nodes(b, quantity, time, level)
                if va.keys() != vb.keys():
                    print(f"t={time} L{level} {quantity}: different nodes")
                    return False
                for key in va:
                    values += 1
                    if va[key] != vb[key]:
                        differ += 1
                        first = first or (time, level, quantity, key, va[key], vb[key])
    print(f"fields: {len(times)} times, {values} values on both levels, {differ} differ")
    if first:
        print(f"  first: t={first[0]} L{first[1]} {first[2]}{first[3]}: {first[4]} vs {first[5]}")
        ok = False

    sets = count_mismatches = 0
    worst = 0.0
    for population in ["protons", "electrons"]:
        with h5py.File(f"{a}/ions_pop_{population}_domain.h5", "r") as f:
            time_keys = sorted(f["t"])
        for time_key in time_keys:
            for level in [0, 1]:
                pa = sorted_particles(a, population, time_key, level)
                pb = sorted_particles(b, population, time_key, level)
                sets += 1
                if pa.shape != pb.shape:
                    count_mismatches += 1
                    print(f"  {population} t={float(time_key):g} L{level}: "
                          f"{pa.shape[0]} vs {pb.shape[0]} particles")
                elif pa.size:
                    worst = max(worst, np.abs(pa - pb).max())
    print(f"particles: {sets} sets, {count_mismatches} count mismatches, max |difference| {worst:.2e}")
    return ok and count_mismatches == 0 and worst <= particle_tolerance


def main():
    argv = sys.argv[1:]
    t2_args = argv[argv.index("--") + 1 :] if "--" in argv else []
    argv = argv[: argv.index("--")] if "--" in argv else argv

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--tiling", action="store_true", help="one patch per level vs tiled")
    mode.add_argument("--ranks", type=int, help="tiled, serial vs this many MPI ranks")
    mode.add_argument("--compare", nargs=2, metavar=("A", "B"), help="two diagnostic directories")
    mode.add_argument("--run", nargs=2, metavar=("PATCH_SIZE", "DIR"), help=argparse.SUPPRESS)
    parser.add_argument("--out", default="phare_outputs/compare_runs", help="diagnostics root")
    parser.add_argument("--particle-tolerance", type=float, default=1e-12)
    cli = parser.parse_args(argv)

    if cli.run:
        run(int(cli.run[0]), cli.run[1], t2_args)
        return 0

    if cli.compare:
        a, b = cli.compare
    elif cli.tiling:
        a, b = f"{cli.out}/one_patch", f"{cli.out}/tiled"
        launch(ONE_PATCH, a, t2_args)
        launch(TILED, b, t2_args)
    else:
        a, b = f"{cli.out}/serial", f"{cli.out}/ranks{cli.ranks}"
        launch(TILED, a, t2_args)
        launch(TILED, b, t2_args, ranks=cli.ranks)

    same = compare(a, b, cli.particle_tolerance)
    print("IDENTICAL" if same else "DIFFERENT")
    return 0 if same else 1


if __name__ == "__main__":
    sys.exit(main())
