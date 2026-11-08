#
# one fixed permute_harris_2d config, repeated, for comparing two builds'
#  PHLOP_SCOPE timers, see compare.sh
#  env PHARE_BENCH_REPEAT=i: only repeat i, compare.sh runs each repeat in its
#   own process: repeats in one process inherit the previous one's memory peak
#   and warm state, later repeats retire fewer instructions
#
"""
PHARE_SCOPE_TIMING=1 PHARE_BENCH_TAR=0 PHARE_BENCH_MON=0 PHARE_BENCH_OUT=.phare_bench/master \\
  python3 -Ou tools/bench/functional/compare_harris_2d.py
python3 -m phlop.timing.scope_compare .phare_bench/master .phare_bench/pr -p 'rank.*.txt'
"""

import os

from pyphare.simulator.simulator import startMPI

repeats = int(os.getenv("PHARE_BENCH_REPEATS", "3"))
repeat = os.getenv("PHARE_BENCH_REPEAT")
steps = int(os.getenv("PHARE_BENCH_STEPS", "10"))
time_step = 0.001

permutables = [
    ("interp_order", [1]),
    ("max_nbr_levels", [3]),
    ("tag_buffer", [3]),
    ("tagging_threshold", [0.1]),
    ("final_time", [steps * time_step]),
    ("repeat", range(repeats) if repeat is None else [int(repeat)]),
]


if __name__ == "__main__":
    from permute_harris_2d import config
    from permutor import execute

    startMPI()
    execute(config, permutables)
