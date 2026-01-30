#
# parsing PHARE scope funtion timers
#

import sys
import argparse
import numpy as np
from dataclasses import dataclass, field

from pyphare.pharesee.run import Run
import phlop.timing.scope_timer as phst

substeps_per_finer_level = 4


@dataclass
class ScopeTimerFile(phst.ScopeTimerFile):
    run: Run
    rank: str
    advances: list = field(default_factory=lambda: [])
    particles_per_level_per_time_step: dict = field(default_factory=lambda: {})

    def __post_init__(self):
        self.advances = [
            root for root in self.roots if self(root.k) == "Simulator::advance"
        ]
        self.pcount_hier, self.particles_per_level_per_time_step = self._particles()

    @property
    def sim(self):
        if sim := self.pcount_hier.sim:
            return sim
        raise ValueError("Simulation is None, check your diagnostics directory")

    def _particles(self):
        """
        Extract particle count per timestep per level from phlop logging
        """
        from pyphare.pharesee.hierarchy.fromh5 import get_times_from_h5

        filepath = self.run.path + "/particle_count.h5"
        all_times = get_times_from_h5(filepath)

        particles_per_level_per_time_step = {}
        pcount_hier = None
        seen_levels = []
        for it, time in enumerate(all_times):
            pcount_hier = self.run.GetParticleCount(time)
            for ilvl, lvl in pcount_hier.levels(time).items():
                if ilvl not in seen_levels:
                    seen_levels += [ilvl]
                pc = sum([p.attrs["particle_count"] for p in lvl.patches])
                if ilvl not in particles_per_level_per_time_step:
                    particles_per_level_per_time_step[ilvl] = np.zeros(
                        len(all_times), dtype=int
                    )
                particles_per_level_per_time_step[ilvl][it] = pc
        return pcount_hier, particles_per_level_per_time_step

    def advance_times_for_L(self, ilvl):
        """
        Extract time in nanoseconds for each substep found for level i
        """
        times = []
        n_levels = self.n_levels()
        assert n_levels == len(self.particles_per_level_per_time_step), "not good"
        for root in self.advances:
            advance_levels = [
                c for c in root.c if self(c.k) == "SolverPPC::advanceLevel"
            ]

            Li_to = 1  # default level 0 value
            for i in range(1, ilvl + 1):
                Li_to += substeps_per_finer_level**i

            Li_from = 0
            for i in range(ilvl):
                Li_from += substeps_per_finer_level**i

            times += [node.t for node in advance_levels[Li_from:Li_to]]

        assert len(times)
        return np.asarray(times, dtype=float)

    def n_levels(self):
        """
        Count number of logs for SolverPPC::advanceLevel per substep per level
        This is mostly to be sure the logs match the number of levels in the diag
        """
        first_advance = self.advances[0]  # all should be the same
        s = sum(1 for c in first_advance.c if self(c.k) == "SolverPPC::advanceLevel")
        if s == 1:
            return 1
        if s == 1 + substeps_per_finer_level:
            return 2
        if s == 1 + substeps_per_finer_level + substeps_per_finer_level**2:
            return 3
        raise ValueError("Add for 4 levels if needed")

    def time_steps_for_L(self, i):
        """
        Get all substep timesteps for the particlular level
        """
        final_time = self.sim.final_time
        time_step = self.sim.time_step
        if i == 0:
            return np.arange(
                time_step, final_time + time_step, final_time / (final_time / time_step)
            )
        sub_steps = self.steps_per_coarse_timestep_for_L(i)
        ts = time_step / sub_steps
        return np.arange(
            ts,
            final_time + ts,
            final_time / (final_time / time_step) / sub_steps,
        )

    def steps_per_coarse_timestep_for_L(self, i):
        return substeps_per_finer_level**i

    def normalised_times_for_L(self, ilvl):
        """
        Normalise substep time against particle count for that level
          at the most recent coarse time, no refined timesteps
        Particle counts may include init dump, so be one bigger.
        """
        times = self.advance_times_for_L(ilvl)
        counts = len(self.particles_per_level_per_time_step[ilvl])

        # trim init particle count for lvl
        Li_times = (
            self.particles_per_level_per_time_step[ilvl]
            if counts == len(times)
            else self.particles_per_level_per_time_step[ilvl][1:]
        )
        if ilvl == 0:
            return times / Li_times
        substeps = self.steps_per_coarse_timestep_for_L(ilvl)
        norm_times = times.copy()
        return (
            norm_times.reshape(int(times.shape[0] / substeps), substeps)
            / Li_times.reshape(Li_times.shape[0], 1)
        ).reshape(times.shape[0])


def file_parser(run, rank, times_filepath):
    supe = phst.file_parser(times_filepath)
    return ScopeTimerFile(supe.id_keys, supe.roots, run, str(rank))


def write_root_as_csv(scope_timer_file, outfile, headers=None, regex=None):
    from contextlib import redirect_stdout

    with open(outfile, "w") as f:
        with redirect_stdout(f):
            print_root_as_csv(scope_timer_file, headers, regex)


def print_root_as_csv(scope_timer_file, n_parts, headers=None, regex=None):
    stf = scope_timer_file  # alias
    stf = file_parser(stf) if isinstance(stf, str) else stf

    if headers:
        print(",".join(headers))
    for root in stf.roots:
        s = stf(root.k)
        if regex and regex not in s:
            continue
        bits = s.split(",")
        print(f"{s}{root.t},{root.t/n_parts}")


def print_variance_across(scope_timer_filepath=None):
    if scope_timer_filepath is None:  # assume cli
        parser = argparse.ArgumentParser()
        parser.add_argument("-f", "--file", default=None, help="timer file")
        scope_timer_filepath = parser.parse_args().file
        if not scope_timer_filepath:
            parser.print_help()
            sys.exit(1)
    phst.print_variance_across(scope_timer_filepath)


def _cli_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", default=None, help="timer file")
    parser.add_argument(
        "-F", "--filter", default=None, help="filter if function supports it"
    )
    return parser


def print_scope_timings(scope_timer_filepath=None, sort_worst_first=True, root_id=None):
    if scope_timer_filepath is None:  # assume cli
        parser = _cli_args()
        args = parser.parse_args()
        scope_timer_filepath = args.file
        if not scope_timer_filepath:
            parser.print_help()
            sys.exit(1)
        if args.filter:
            root_id = args.filter
    phst.print_scope_timings(scope_timer_filepath, sort_worst_first, root_id)


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
