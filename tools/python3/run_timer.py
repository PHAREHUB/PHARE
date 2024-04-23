#
# parsing PHARE scope funtion timers
#

import numpy as np
from dataclasses import dataclass, field

from pyphare.pharesee.run import Run
from pyphare.pharesee.hierarchy import hierarchy_from


@dataclass
class RunTimerNode:
    k: int
    t: int
    c: list = field(default_factory=lambda: [])


@dataclass
class RunTimerFile:
    run: Run
    rank: str
    id_keys: dict = field(default_factory=lambda: {})
    roots: list = field(default_factory=lambda: [])
    advances: list = field(default_factory=lambda: [])
    particles_per_level_per_time_step: dict = field(default_factory=lambda: {})

    def fn_for(self, id):
        return self.id_keys[id]

    def __call__(self, id):
        return self.fn_for(id)

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
    # compressed form to reduce string usage

    id_keys = {}
    roots = []
    curr = None
    c_depth = 0
    stack = {i: 0 for i in range(128)}
    stack_size = 0

    def _parent():
        assert stack_size >= 0
        n = roots[stack[0]]
        for s in range(1, stack_size):
            assert stack[s] >= 0
            n = n.c[stack[s]]
        return n

    with open(times_filepath, "r") as file:
        while line := file.readline():
            line = line.rstrip()
            if not line:
                break
            bits = line.split(" ")
            id_keys[bits[0]] = bits[1]
        while line := file.readline():
            line = line.rstrip()  # drop new line characters
            stripped_line = line.strip()
            if not stripped_line:  # last line might be blank
                continue
            idx = len(line) - len(stripped_line)  # how many space indents from left
            bits = stripped_line.split(" ")
            node = RunTimerNode(*bits)
            if idx == 0:  # is root node
                stack[0] = len(roots)
                roots.append(node)
                stack_size = 0
                c_depth = 0
            else:  # is not root node
                if idx > c_depth:
                    parent = curr
                    curr.c.append(node)
                    stack_size += 1
                elif idx == c_depth:
                    parent = _parent()
                    parent.c.append(node)
                elif idx < c_depth:
                    stack_size -= c_depth - idx
                    parent = _parent()
                    parent.c.append(node)
                stack[stack_size] = len(parent.c) - 1
                c_depth = idx
            curr = node
    return RunTimerFile(run, str(rank), id_keys, roots)
