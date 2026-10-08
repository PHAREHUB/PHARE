#
#
#
import os
import sys
import json
import shutil
import datetime
import itertools
import subprocess
from pathlib import Path

from pyphare import cpp
import pyphare.pharein as ph
from pyphare.simulator.monitoring import MonitoringOptions
from pyphare.simulator.simulator import Simulator


def datetime_now():
    return datetime.datetime.now().replace(microsecond=0).strftime("%Y%m%d-%H%M%S")


MAKE_TAR_FILE = True
log_dir = Path(".log")
local_dir = Path(".phare")
out_dir = Path(f".phare_bench/{datetime_now()}")
monitoring_options = MonitoringOptions(interval=5, rank_modulo=1)


def sort_summaries(path):
    import re

    with open(path) as f:
        lines = [line.rstrip("\n") for line in f if line.strip()]

    def key(line):
        tb = float(re.search(r"'tag_buffer': (\d+)", line).group(1))
        tt = float(re.search(r"'tagging_threshold': ([0-9.]+)", line).group(1))
        ts = float(re.search(r"'tile_size': ([0-9.]+)", line).group(1))
        return (tb, tt, ts)

    lines.sort(key=key)

    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def mkdir(dir):
    dir.mkdir(parents=True, exist_ok=True)


def clean_dir(dir):
    if dir.exists():
        shutil.rmtree(str(dir))


def make_tarfile(outdir, outfile=None, rank0=None):
    if MAKE_TAR_FILE and (cpp.mpi_rank() == 0 if rank0 is None else rank0):
        path = outdir.parent
        outfile = outfile or outdir.name
        args = ["tar", "czf", f"{path/outfile}.tar.gz", "-C", str(path), outdir.name]
        ec = subprocess.call(args)
        assert ec == 0, f"{ec}"
        clean_dir(outdir)


def post_sim(summary, kwargs):
    log = local_dir / "summary.txt"
    with open(log, "w") as f:
        f.write(f"{kwargs}\n")
        f.write(summary)
    to = out_dir / datetime_now()
    shutil.copytree(log_dir, to / log_dir.name)
    shutil.copytree(local_dir, to / local_dir.name)
    make_tarfile(to)


def execute_permutation(config_fn, keys, permutation):
    from pyphare import cpp

    def run():
        ph.global_vars.sim = None
        kwargs = {keys[i]: permutation[i] for i in range(len(keys))}
        sim = Simulator(config_fn(**kwargs)).run(monitoring=monitoring_options)
        summary = sim.summary()
        sim.reset()
        cpp.mpi_barrier()
        if cpp.mpi_rank() == 0:
            post_sim(summary, kwargs)

    run()
    cpp.mpi_barrier()


def execute(config_fn, permutables):
    from pyphare import cpp

    if cpp.mpi_rank() == 0:
        mkdir(out_dir)  # should be unique
        mkdir(log_dir)
        mkdir(local_dir / "timings")
    cpp.mpi_barrier()
    keys = tuple(e[0] for e in permutables)
    for permutation in itertools.product(*(e[1] for e in permutables)):
        execute_permutation(config_fn, keys, permutation)
    make_tarfile(out_dir)


# Each permutation runs in its own `mpirun -n N python3 <script>` subprocess, so MPI ranks and
# thread pools (PHARE_THREAD_POOLS * PHARE_THREADS_PER_POOL, read once per process) can be varied
# per run. `launches` are dicts of {"mpirun": ranks, "pools": n_pools, "threads": per_pool}.
# Permutable keys starting with "PHARE_" are env vars for the subprocess (e.g.
# "PHARE_TILING_MAX_TILE_SIZE"), not config kwargs. skip(record) -> True drops a combination.

CHILD_ARG = "--permutation"


def _split_env(kwargs):
    env = {k: v for k, v in kwargs.items() if k.startswith("PHARE_")}
    return {k: v for k, v in kwargs.items() if k not in env}, env


def execute_launches(config_fn, permutables, launches, skip=None, mpirun="mpirun"):
    if CHILD_ARG in sys.argv:
        return _execute_child(config_fn)

    mkdir(out_dir)
    keys = tuple(e[0] for e in permutables)
    script = Path(sys.argv[0]).resolve()

    for permutation in itertools.product(*(e[1] for e in permutables)):
        for launch in launches:
            record = {**dict(zip(keys, permutation)), **launch}
            if skip and skip(record):
                continue
            _, env_kwargs = _split_env(record)
            env = {
                "PHARE_SCOPE_TIMING": "1",  # timers are what we're here for, env can override
                **os.environ,
                **{k: str(v) for k, v in env_kwargs.items()},
                "PHARE_THREAD_POOLS": str(launch["pools"]),
                "PHARE_THREADS_PER_POOL": str(launch["threads"]),
            }
            cmd = [mpirun, "-n", str(launch["mpirun"]), sys.executable, "-u", str(script)]
            cmd += [CHILD_ARG, json.dumps(record), "--out-dir", str(out_dir)]
            print(" ".join(cmd[:4]), record, flush=True)
            ec = subprocess.call(cmd, env=env)
            if ec != 0:
                print(f"FAILED ({ec}): {record}", flush=True)

    make_tarfile(out_dir, rank0=True)  # no MPI in this (parent) process


def _execute_child(config_fn):
    global out_dir
    from pyphare import cpp
    from pyphare.simulator.simulator import startMPI

    startMPI()
    record = json.loads(sys.argv[sys.argv.index(CHILD_ARG) + 1])
    out_dir = Path(sys.argv[sys.argv.index("--out-dir") + 1])

    if cpp.mpi_rank() == 0:  # only this run's timers/stats may end up in its tarball
        clean_dir(log_dir)
        clean_dir(local_dir)
        mkdir(log_dir)
        mkdir(local_dir / "timings")
    cpp.mpi_barrier()

    launch_keys = ("mpirun", "pools", "threads")
    kwargs, _ = _split_env({k: v for k, v in record.items() if k not in launch_keys})

    ph.global_vars.sim = None
    sim = Simulator(config_fn(**kwargs)).run(monitoring=monitoring_options)
    summary = sim.summary()
    sim.reset()
    cpp.mpi_barrier()
    if cpp.mpi_rank() == 0:
        post_sim(summary, record)
    cpp.mpi_barrier()
