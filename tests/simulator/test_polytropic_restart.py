#!/usr/bin/env python3
#
# restart of a run with the polytropic electron closure and a refined box. A bump of Pe
# is convected by a uniform flow, so that the pressure carries a state that a restart must
# bring back: a restart that seeded the entropy again from the user's profiles would put
# the bump back where it started. The run restarted from t = final_time / 2 must give, at
# final_time, the same pressure as the uninterrupted one on every level, and the bump
# must have kept on moving.

import copy
import unittest

import numpy as np

import pyphare.pharein as ph
from pyphare import cpp
from pyphare.pharesee.run import Run
from pyphare.simulator.simulator import Simulator, startMPI

ph.NO_GUI()

cells = 100
dl = 0.5
V = 2.0
P0 = 0.5
amplitude = 0.5
width = 2.5
x0 = 8.0
gamma = 5.0 / 3.0
time_step = 0.01
final_time = 5.0
restart_time = 2.5
out = "phare_outputs/polytropic_restart/"


def Pe(x):
    return P0 * (1.0 + amplitude * np.exp(-(((x - x0) / width) ** 2)))


def config(diag_dir, restart_options):
    sim = ph.Simulation(
        time_step=time_step,
        final_time=final_time,
        cells=cells,
        dl=dl,
        gamma=gamma,
        largest_patch_size=20,
        smallest_patch_size=20,
        refinement_boxes={"L0": {"B0": [(30,), (70,)]}},
        diag_options={
            "format": "phareh5",
            "options": {"dir": diag_dir, "mode": "overwrite"},
        },
        restart_options=restart_options,
    )

    def one(x):
        return 1.0

    def zero(x):
        return 0.0

    def vth(x):
        return 0.3

    ph.MaxwellianFluidModel(
        bx=one,
        by=zero,
        bz=zero,
        protons={
            "charge": 1,
            "density": one,
            "vbulkx": lambda x: V,
            "vbulky": zero,
            "vbulkz": zero,
            "vthx": vth,
            "vthy": vth,
            "vthz": vth,
            "nbr_part_per_cell": 200,
            "init": {"seed": 2026},
        },
    )
    ph.ElectronModel(closure="polytropic", Pe=Pe)
    stamps = [restart_time, final_time]
    ph.FluidDiagnostics(quantity="electron_pressure", write_timestamps=stamps)
    ph.FluidDiagnostics(quantity="charge_density", write_timestamps=stamps)
    return sim


def domain_nodes(hier, time, ilvl):
    xs, vs = [], []
    for patch in hier.level(ilvl, time).patches:
        pdata = next(iter(patch.patch_datas.values()))
        g = pdata.ghosts_nbr[0]
        domain = slice(g, len(pdata.x) - g)
        xs.append(np.asarray(pdata.x)[domain])
        vs.append(np.asarray(pdata.dataset[:])[domain])
    x, idx = np.unique(np.round(np.concatenate(xs), 9), return_index=True)
    return x, np.concatenate(vs)[idx]


class PolytropicRestartTest(unittest.TestCase):
    def tearDown(self):
        ph.global_vars.sim = None

    def test_restart(self):
        dir0, dir1 = out + "continuous", out + "restarted"

        restart_options = {
            "dir": dir0,
            "mode": "overwrite",
            "timestamps": [restart_time],
        }
        Simulator(config(dir0, restart_options)).run().reset()
        ph.global_vars.sim = None

        restarted = copy.deepcopy(restart_options)
        restarted.update(restart_time=restart_time)
        del restarted["timestamps"]
        Simulator(config(dir1, restarted)).run().reset()

        run0, run1 = Run(dir0), Run(dir1)
        for ilvl in (0, 1):
            x0_, pe0 = domain_nodes(
                run0.GetPe(final_time, all_primal=False), final_time, ilvl
            )
            x1_, pe1 = domain_nodes(
                run1.GetPe(final_time, all_primal=False), final_time, ilvl
            )
            np.testing.assert_array_equal(x0_, x1_)
            # exactly the same in a serial run; with several ranks the sums of the
            # deposited moments are not made in the same order after a restart, which
            # changes the last bit of the float32 diagnostics
            rtol = 0 if cpp.mpi_size() == 1 else 1e-6
            np.testing.assert_allclose(
                pe1,
                pe0,
                rtol=rtol,
                atol=0,
                err_msg=f"level {ilvl} differs after a restart",
            )

        # the bump went on moving: at V t from where it started, not back at the start
        x, pe = domain_nodes(run1.GetPe(final_time, all_primal=False), final_time, 0)
        _, n = domain_nodes(run1.GetNi(final_time), final_time, 0)
        peak = x[
            np.argmax(pe / n**gamma)
        ]  # the entropy: Pe has the noise of the density
        print(f"bump peak at {peak}, expected about {x0 + V * final_time}")
        self.assertLess(abs(peak - (x0 + V * final_time)), 3.0)


if __name__ == "__main__":
    startMPI()
    unittest.main()
