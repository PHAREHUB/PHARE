import yaml
from pathlib import Path


def get_yaml_for(file):
    with open(file) as stream:
        try:
            return yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)


files = list(Path(".phare/stats").glob("*.yaml"))
roots = [get_yaml_for(file) for file in files]
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
