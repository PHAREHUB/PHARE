#
# parsing PHARE scope funtion timers
#

import sys
import argparse
import numpy as np
from dataclasses import dataclass, field


@dataclass
class ScorePScoreUserLine:
    type: str
    max_buf: float
    visits: float
    time_second: float
    time_percent: float
    time_per_visit: float
    region: str

    def __lt__(self, that):
        return self.time_percent > that.time_percent

    def __repr__(self):
        return f"% {self.time_percent} in: {' '.join(self.region)}"

    def __str__(self):
        return self.__repr__()

    @staticmethod
    def FROM_LINE(line):
        bits = [bit for bit in line.split(" ") if bit]
        return ScorePScoreUserLine(
            bits[0],
            float(bits[1].replace(",", "")),
            float(bits[2].replace(",", "")),
            float(bits[3].replace(",", "")),
            float(bits[4]),
            float(bits[5]),
            bits[6:],
        )


@dataclass
class ScorePScore:
    header: str
    summary: str
    usr_lines: list

    def __post_init__(self):
        ...

    def parsed_user_lines(self, worst_first):
        lines = []
        for line in self.usr_lines:
            lines.append(ScorePScoreUserLine.FROM_LINE(line))

        if worst_first:
            return sorted(lines)

    def print(self, worst_first=True):
        print("".join(self.header) + "\n")
        print("".join(self.summary) + "\n")

        usr_lines = self.parsed_user_lines(worst_first)
        print("".join(str(line) for line in usr_lines))


def scorep_score(filepath):
    section = 0

    lines = [[], [], []]

    with open(filepath) as f:
        line = f.readline().strip()
        if line:
            lines[0].append(line)
        while line := f.readline():
            if not line.strip():
                section += 1
            else:
                lines[section].append(line)

    return ScorePScore(*lines)


def print_score_p_score(filepath=None, worst_first=True):
    if filepath is None:  # assume cli
        parser = _cli_args()
        args = parser.parse_args()
        filepath = args.file
        if not filepath:
            parser.print_help()
            sys.exit(1)
    scorep_score(filepath).print(worst_first)


def _cli_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", default=None, help="timer file")
    return parser


if __name__ == "__main__":
    if len(sys.argv) == 1:
        print("usage: $function_name -h")
        print(
            "available functions:\n\t"
            + "\n\t".join([k for k, v in globals().items() if k.startswith("print_")]),
        )
    elif len(sys.argv) > 1:
        fn = sys.argv[1]
        sys.argv = [sys.argv[0]] + sys.argv[2:]
        if fn not in globals():
            raise ValueError("requested function does not exist")
        globals()[fn]()
