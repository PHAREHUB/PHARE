#
# parsing PHARE scope funtion timers
#

import os
import sys
import json
import tarfile
import argparse
import numpy as np
from io import BytesIO
from pathlib import Path
from typing import Callable
from dataclasses import asdict, dataclass, field

from pyphare.pharesee.run import Run
import phlop.timing.scope_timer as phst


def _as_lines(input):
    if type(input) is bytes:
        return input.decode("utf8").split("\n")
    with open(input, "r") as file:
        return file.readlines()


def cores(n=1, c=1):
    return n * c


def steps_per(ilvl):
    return 4**ilvl


@dataclass
class LevelStats:  # one level of Simulator.summary(), see model_stats.hpp
    level: int
    model: str  # "hybrid" or "mhd"
    patches: int
    cells: int
    patch_ghosts: int  # ghost cells covered by a same level neighbour
    level_ghosts: int  # the other ghost cells
    remote_ghosts: int  # the patch ghosts owned by another rank
    remote_ranks_max: int  # other ranks owning patch ghosts, max over ranks
    remote_ranks_avg: float
    cell_imbalance: float  # max per rank over average, 1 is perfectly balanced
    particles: int | None = None  # hybrid levels only
    particle_imbalance: float | None = None  # hybrid levels only

    @property
    def is_hybrid(self):
        return self.model == "hybrid"

    @property
    def cell_ratio(self):  # ghost over domain cells
        return (self.patch_ghosts + self.level_ghosts) / self.cells

    @property
    def remote_pc(self):  # % of patch ghosts owned by another rank
        return 100 * self.remote_ghosts / self.patch_ghosts if self.patch_ghosts else 0


@dataclass
class SummaryFile:
    kwargs: dict
    samples: list  # [(time, [LevelStats, ...]), ...] in time order

    @property
    def levels(self):  # the final hierarchy
        return self.samples[-1][1]

    @property
    def particles(self):  # None without hybrid levels
        hybrid = [lvl.particles for lvl in self.levels if lvl.is_hybrid]
        return sum(hybrid) if hybrid else None

    def normalized_nanoseconds_for_time_per_push(self, time_ns, n_cores):
        """None without particles, as for MHD only runs"""
        ppp = sum(
            lvl.particles * steps_per(lvl.level) for lvl in self.levels if lvl.is_hybrid
        )
        if ppp == 0:
            return None
        return ppp, time_ns / ppp * n_cores


def parse_summary_file(input):
    """the permutor's summary.json: the run's kwargs and Simulator.summary() samples"""
    summary = json.loads(input if type(input) is bytes else Path(input).read_text())
    return SummaryFile(
        summary["kwargs"],
        [
            (sample["time"], [LevelStats(**lvl) for lvl in sample["levels"]])
            for sample in summary["samples"]
        ],
    )


def percent_level_refined(levels, ndim):
    """per level, % of it covered by the next finer level, None on the finest"""
    return [
        (
            int(levels[i + 1].cells / 2**ndim / levels[i].cells * 100)
            if i + 1 < len(levels)
            else None
        )
        for i in range(len(levels))
    ]


def average_L0_advance_time_from(timer_files):
    advance_times = [tf.advances[0].run_time for tf in timer_files]
    return np.average(advance_times)


def max_memory_required_from(stat_files):
    ret = 0
    for file in stat_files:
        for snapshot in file.snapshots:
            ret = max(ret, snapshot.mem_used_mb)
    return ret


def _print_stats(summary_file, timer_files, stat_files):
    """one JSON object per run on stdout, read back with parse_stats_line"""
    if not timer_files:
        print(
            "no timer files, run with PHARE_SCOPE_TIMING=1",
            summary_file.kwargs,
            file=sys.stderr,
        )
        return
    time = average_L0_advance_time_from(timer_files)
    pushes = summary_file.normalized_nanoseconds_for_time_per_push(
        time, len(timer_files)
    )
    ndim = len(summary_file.kwargs["cells"])
    samples = [
        dict(time=t, levels=_level_records(levels, ndim))
        for t, levels in summary_file.samples
    ]
    record = dict(
        kwargs=summary_file.kwargs,
        time_s=time / 1e9,  # average L0 advance
        mem_mb=max_memory_required_from(stat_files),
        particles=summary_file.particles,
        pushes=pushes[0] if pushes else None,
        ppp_ns=pushes[1] if pushes else None,  # per particle push, normalised by cores
        levels=samples[-1]["levels"],  # the final hierarchy
        samples=samples,
    )
    print(json.dumps(record))


def parse_stats_line(line):
    """a print_summary output line as a dict, None for anything else"""
    try:
        record = json.loads(line)
    except json.JSONDecodeError:
        return None
    return record if isinstance(record, dict) and "levels" in record else None


def extract_nested_tarfile(data):
    with tarfile.open(fileobj=data, mode="r:*") as nestedTar:

        def extract(member):
            if file := nestedTar.extractfile(member):
                return file.read()

        return {
            k: v
            for k, v in {member.name: extract(member) for member in nestedTar}.items()
            if v
        }


def handle_timer_file(data):
    import tools.python3.phloping as phloping

    return phloping.make_scope_timer_file_from(_as_lines(data))


def handle_stat_file(data):
    import tools.python3.stats as stats

    return stats.lines_parser(data)


def handle_pertumation(objects):
    stat_files = [
        handle_stat_file(v) for k, v in objects.items() if "phare/stats/" in k
    ]
    timer_files = [
        handle_timer_file(v) for k, v in objects.items() if "phare/timings/" in k
    ]
    for k, v in objects.items():
        if k.endswith("phare/summary.json"):
            _print_stats(parse_summary_file(v), timer_files, stat_files)


def extract_nested_directory(path):
    def extract(file):
        with open(file, "rb") as f:
            return f.read()

    return {
        str(Path(d) / file): extract(str(Path(d) / file))
        for d, _, files in os.walk(path)  # Path.walk needs python 3.12
        for file in files
    }


def handle_nested_directory(path):
    handle_pertumation(extract_nested_directory(path))


def handle_nested_tarfile(nested):
    if type(nested) is not BytesIO:
        with open(str(nested), "rb") as f:
            nested = BytesIO(f.read())
    handle_pertumation(extract_nested_tarfile(nested))


def handle_tarfile(path):
    with tarfile.open(str(path), "r:gz") as tar:
        for member in tar:
            if member.name.endswith(".tar.gz"):
                handle_nested_tarfile(BytesIO(tar.extractfile(member).read()))


def this_or_that(pred, a, b, *args, **kwargs):
    return a(*args, **kwargs) if pred else b(*args, **kwargs)


def is_tar_file(any):
    return str(any).endswith(".tar.gz")


def handle_directory(path):
    path = Path(path)
    if (path / "summary.json").exists():  # a single run's .phare
        return handle_nested_directory(path)
    for d in path.iterdir():  # e.g. .phare_bench: outer archives or run directories
        this_or_that(is_tar_file(d), handle_tarfile, handle_nested_directory, d)


def handle_input(path):
    this_or_that(is_tar_file(path), handle_tarfile, handle_directory, path)


def print_summary(filepath=None, root_id=None):
    if filepath:
        return handle_input(filepath)

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", default=None, help="input")
    args = parser.parse_args()
    if args.input:
        return handle_input(args.input)
    else:
        parser.print_help()


def _level_records(levels, ndim):  # LevelStats as dicts, with the computed values
    refined = percent_level_refined(levels, ndim)
    return [
        dict(
            asdict(lvl),
            cell_ratio=lvl.cell_ratio,
            remote_pc=lvl.remote_pc,
            refined_pc=refined[i],
        )
        for i, lvl in enumerate(levels)
    ]


if __name__ == "__main__":
    if len(sys.argv) == 1:
        print("usage: $function_name -h")
        print(
            "available functions:\n\t"
            + "\n\t".join([k for k, v in globals().items() if k.startswith("print")]),
        )
    elif len(sys.argv) > 1:
        fn = sys.argv[1]
        sys.argv = [sys.argv[0]] + sys.argv[2:]
        if fn not in globals():
            raise ValueError("requested function does not exist")
        try:
            globals()[fn]()
        except BrokenPipeError:
            ...
