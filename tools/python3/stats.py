import yaml
from typing import List
from pathlib import Path
from dataclasses import dataclass, field


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
    snapshots: List[StatFileSnapshot]
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
    roots = [get_yaml_for(file) for file in files]

    if not roots:
        raise ValueError("No snapshots found!")

    # ranks sample independently, counts can differ: sum what each sample index has
    snapshots = [root.get("snapshots") or [] for root in roots]
    times = max(len(s) for s in snapshots)
    if not times:
        print("mem: unavailable, no snapshots")
        return

    mem = [0] * times

    for snaps in snapshots:
        for tidx, snapshot in enumerate(snaps):
            bits = snapshot["v"].split(",")
            mem[tidx] += int(bits[1])

    print("mem: ", mem, " max(", max(mem), ")")


if __name__ == "__main__":
    main()
