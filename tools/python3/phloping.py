#
# parsing PHARE scope funtion timers
#

import numpy as np
from dataclasses import dataclass, field

from pyphare.pharesee.run import Run
from pyphare.pharesee.hierarchy import hierarchy_from

from phlop.timing.scope_timer import ScopeTimerFile as phScopeTimerFile
from phlop.timing.scope_timer import file_parser as phfile_parser


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
        self.particles_per_level_per_time_step = self._particles()

    @property
    def sim(self):
        if sim := self.run.GetAllAvailableQties().sim:
            return sim
        raise ValueError("Simulation is None, check your diagnostics directory")

    def _particles(self):
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
        return particles_per_level_per_time_step

    def advance_times_for_L(self, i):
        times = []
        n_levels = self.n_levels()
        assert n_levels == len(self.particles_per_level_per_time_step), "not good"
        for root in self.advances:
            advance_levels = [
                c for c in root.c if self(c.k) == "SolverPPC::advanceLevel"
            ]
            if i == 0:
                times += [advance_levels[0].t]
            elif i == 1 and n_levels == 2:
                times += [node.t for node in advance_levels[1:5]]
            elif i == 1 and n_levels == 3:
                raise RuntimeError("not finished")
        assert len(times)
        return np.asarray(times, dtype=float)

    def n_levels(self):
        _0 = self.advances[0]  # all should be the same
        s = sum(1 for c in _0.c if self(c.k) == "SolverPPC::advanceLevel")
        substeps_per_finer_level = 4
        if s == 1:
            return 1
        if s == 1 + substeps_per_finer_level:
            return 2
        if s == 1 + substeps_per_finer_level + substeps_per_finer_level**2:
            return 3
        raise RuntimeError("finish for L4+")

    def time_steps_for_L(self, i):
        final_time = self.sim.final_time
        time_step = self.sim.time_step
        if i == 0:
            return np.arange(
                time_step, final_time + time_step, final_time / (final_time / time_step)
            )
        substeps_per_finer_level = 4
        ts = time_step / substeps_per_finer_level**i
        return np.arange(
            ts,
            final_time + ts,
            final_time / (final_time / time_step) / substeps_per_finer_level,
        )

    def normalised_times_for_L(self, ilvl):
        times = self.advance_times_for_L(ilvl)
        if ilvl == 0:
            return times / self.particles_per_level_per_time_step[0]
        elif ilvl == 1:
            norm_times = times.copy()
            return (
                norm_times.reshape(int(times.shape[0] / 4), 4)
                / self.particles_per_level_per_time_step[1].reshape(
                    self.particles_per_level_per_time_step[1].shape[0], 1
                )
            ).reshape(times.shape[0])
        else:
            raise RuntimeError("not finished")


def file_parser(run, rank, times_filepath):
    supe = phfile_parser(times_filepath)
    return ScopeTimerFile(supe.id_keys, supe.roots, run, str(rank))
