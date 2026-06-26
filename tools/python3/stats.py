import yaml
from typing import List
from pathlib import Path
from dataclasses import dataclass, field


@dataclass
class StatFileSnapshot:
    open_files: int
    mem_used_mb: int
    cpu_usage: float

    @staticmethod
    def FROM_CSV(input):
        a, b, c = input.split(",")
        return StatFileSnapshot(int(a), int(b), float(c))


@dataclass
class StatFile:
    rank: int
    cli_args: map
    headers: list
    start: str
    end: str
    snapshots: List[StatFileSnapshot]

    def __post_init__(self):
        if self.snapshots and type(self.snapshots[0]) is not StatFileSnapshot:
            self.snapshots = [StatFileSnapshot.FROM_CSV(v["v"]) for v in self.snapshots]


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

    times = len(roots[0]["snapshots"])
    consistent = all(times == len(root["snapshots"]) for root in roots)

    if not consistent:
        raise ValueError("Different number of snapshots found across files!")

    mem = [0] * times
    cpu = [0] * times

    for root in roots:
        for tidx, snapshot in enumerate(root["snapshots"]):
            bits = snapshot["v"].split(",")
            mem[tidx] += int(bits[1])

    print("mem: ", mem, " max(", max(mem), ")")


if __name__ == "__main__":
    main()
