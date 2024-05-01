#
# parsing PHARE scope funtion timers
#

import numpy as np
from dataclasses import dataclass, field

from pyphare.pharesee.run import Run
from pyphare.pharesee.hierarchy import hierarchy_from

from phlop.timing.scope_timer import ScopeTimerFile as phScopeTimerFile
from phlop.timing.scope_timer import file_parser as phfile_parser


substeps_per_finer_level = 4


@dataclass
class ScopeTimerFile(phScopeTimerFile):
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
        bad_diag_error = (
            "Simulation was not configured with Particle Count info diagnostic"
        )
        pcount_hier = hierarchy_from(h5_filename=self.run.path + "/particle_count.h5")
        particles_per_level_per_time_step = {  # per coarse timestep only
            li: np.zeros(len(self.advances))
            for li in range(len(pcount_hier.data_files["t"][pcount_hier.times()[0]]))
        }
        for ti, t in enumerate(pcount_hier.times()[1:]):
            for plk, pl in pcount_hier.data_files["t"][t].items():
                pc = 0
                for pid, p in pl.items():
                    if "particle_count" not in p.attrs:
                        raise ValueError(bad_diag_error)
                    if pid.split("#")[0][1:] == self.rank:
                        pc += p.attrs["particle_count"]
                if pc == 0:
                    pc += 1  # avoid div by 0 for rank missing patch on level
                particles_per_level_per_time_step[int(plk[2:])][ti] += pc
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
        """
        times = self.advance_times_for_L(ilvl)
        if ilvl == 0:
            return times / self.particles_per_level_per_time_step[0]
        substeps = self.steps_per_coarse_timestep_for_L(ilvl)
        norm_times = times.copy()
        return (
            norm_times.reshape(int(times.shape[0] / substeps), substeps)
            / self.particles_per_level_per_time_step[ilvl].reshape(
                self.particles_per_level_per_time_step[ilvl].shape[0], 1
            )
        ).reshape(times.shape[0])


def file_parser(run, rank, times_filepath):
    supe = phfile_parser(times_filepath)
    return ScopeTimerFile(supe.id_keys, supe.roots, run, str(rank))
