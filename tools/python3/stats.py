from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import yaml


@dataclass
class StatFileSnapshot:
    open_files: int
    mem_used_mb: int
    cpu_usage: float
    extra: dict = field(default_factory=dict)  # other columns by header, e.g. gpu0_util

    @staticmethod
    def FROM_CSV(input, headers=("fds", "mem_used_mb", "cpu_usage")):
        values = dict(zip(headers, input.split(",")))
        return StatFileSnapshot(
            int(values.pop("fds")),
            int(values.pop("mem_used_mb")),
            float(values.pop("cpu_usage")),
            values,
        )


@dataclass
class StatFile:
    rank: int
    cli_args: map
    headers: list
    start: str
    end: str
    snapshots: list[StatFileSnapshot]
    gpus: list = field(default_factory=list)  # phlop with GPUs found

    def __post_init__(self):
        # "snapshots:" without value loads as None
        self.snapshots = self.snapshots or []
        if self.snapshots and type(self.snapshots[0]) is not StatFileSnapshot:
            self.snapshots = [
                StatFileSnapshot.FROM_CSV(v["v"], self.headers) for v in self.snapshots
            ]


def lines_parser(input):
    return StatFile(**yaml.safe_load(input))


def get_yaml_for(file):
    with open(file) as stream:
        try:
            return yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            raise RuntimeError(f"Failed to parse {file}") from exc


def main():  # default
    files = list(Path(".phare/stats").glob("*.yaml"))
    if not files:
        raise ValueError("No snapshots found!")

    stat_files = [StatFile(**get_yaml_for(file)) for file in files]
    sampled = [sf for sf in stat_files if sf.snapshots]
    if not sampled:
        print("mem: unavailable, no snapshots")  # monitoring disabled
        return

    # ranks sample on their own clocks, at their own intervals: put every rank on one time
    # grid, each with its reading nearest each grid time, so the total sums values from about
    # the same time
    times = [_sample_times(sf) for sf in sampled]
    intervals = sorted({sf.cli_args["interval"] for sf in sampled})
    first, last = min(t[0] for t in times), max(t[-1] for t in times)
    # up to the last sample
    grid = np.arange(first, last + intervals[0] / 2, intervals[0])
    total = np.zeros(len(grid))
    peaks = np.empty(len(sampled))
    for i, (sf, t) in enumerate(zip(sampled, times)):
        mb = np.array([s.mem_used_mb for s in sf.snapshots])
        # latest sample up to half a step after the grid time: start times differ by microseconds
        near = np.searchsorted(t, grid + intervals[0] / 2, side="right") - 1
        total += mb[np.clip(near, 0, len(mb) - 1)]
        peaks[i] = mb.max()

    peak = total.argmax()
    print(f"mem over {len(sampled)} sampled rank(s), interval(s) {intervals}s:")
    print(
        f"  job peak ~{total[peak]:.0f} MB at +{grid[peak] - grid[0]:.0f}s"
        f" (rank readings up to {intervals[-1]}s apart)"
    )
    print(
        f"  rank peaks (MB): min {peaks.min():.0f}  median {np.median(peaks):.0f}"
        f"  max {peaks.max():.0f} (rank {sampled[peaks.argmax()].rank})"
    )


def _sample_times(sf):
    """phlop samples min(0.1, interval) after start, then every interval"""
    interval = sf.cli_args["interval"]
    start = sf.start.timestamp()  # yaml loads the ISO time as a datetime
    return start + min(0.1, interval) + interval * np.arange(len(sf.snapshots))


if __name__ == "__main__":
    main()
