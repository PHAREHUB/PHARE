#!/usr/bin/env python3

# basically a harris run

import unittest
from ddt import data, ddt, unpack
import pyphare.pharein as ph
from pyphare.pharesee.hierarchy.hierarchy_utils import hierarchy_compare
from pyphare.pharesee.particles import single_patch_per_level_per_pop_from

from pyphare.simulator.simulator import Simulator, startMPI
from tests.simulator import SimulatorTest

import numpy as np
import matplotlib as mpl

mpl.use("Agg")

from pyphare.cpp import cpp_lib

cpp = cpp_lib()
startMPI()

ndim = 2
interp = 1
mpi_size = cpp.mpi_size()
time_step_nbr = 3
time_step = 0.001
cells = (100, 100)
dl = (0.2, 0.2)
diag_outputs = "phare_outputs/load_balancing_2d"
timestamps = [x * time_step for x in range(time_step_nbr + 1)]


def config(diag_dir, loadbalancing={}):
    sim = ph.Simulation(
        time_step_nbr=time_step_nbr,
        time_step=time_step,
        cells=cells,
        dl=dl,
        refinement="tagging",
        max_nbr_levels=2,
        hyper_resistivity=0.001,
        resistivity=0.001,
        diag_options={
            "format": "phareh5",
            "options": {"dir": diag_dir, "mode": "overwrite"},
        },
    )

    def density(x, y):
        L = sim.simulation_domain()[1]
        return (
            0.2
            + 1.0 / np.cosh((y - L * 0.3) / 0.5) ** 2
            + 1.0 / np.cosh((y - L * 0.7) / 0.5) ** 2
        )

    def densityHigh(x, y):
        rho = y.copy()
        rho[:] = 1
        return rho

    def densityLow(x, y):
        assert cells == (100, 100)  # next lines are dependent
        rho = y.copy()
        rho[:] = 0.1
        rho[np.where(np.isclose(y, 10, atol=0.1))] = 0
        return rho

    def S(y, y0, l):
        return 0.5 * (1.0 + np.tanh((y - y0) / l))

    def by(x, y):
        Lx = sim.simulation_domain()[0]
        Ly = sim.simulation_domain()[1]
        w1 = 0.2
        w2 = 1.0
        x0 = x - 0.5 * Lx
        y1 = y - 0.3 * Ly
        y2 = y - 0.7 * Ly
        w3 = np.exp(-(x0 * x0 + y1 * y1) / (w2 * w2))
        w4 = np.exp(-(x0 * x0 + y2 * y2) / (w2 * w2))
        w5 = 2.0 * w1 / w2
        return (w5 * x0 * w3) + (-w5 * x0 * w4)

    def bx(x, y):
        Lx = sim.simulation_domain()[0]
        Ly = sim.simulation_domain()[1]
        w1 = 0.2
        w2 = 1.0
        x0 = x - 0.5 * Lx
        y1 = y - 0.3 * Ly
        y2 = y - 0.7 * Ly
        w3 = np.exp(-(x0 * x0 + y1 * y1) / (w2 * w2))
        w4 = np.exp(-(x0 * x0 + y2 * y2) / (w2 * w2))
        w5 = 2.0 * w1 / w2
        v1 = -1
        v2 = 1.0
        return (
            v1
            + (v2 - v1) * (S(y, Ly * 0.3, 0.5) - S(y, Ly * 0.7, 0.5))
            + (-w5 * y1 * w3)
            + (+w5 * y2 * w4)
        )

    def bz(x, y):
        return 0.0

    def b2(x, y):
        return bx(x, y) ** 2 + by(x, y) ** 2 + bz(x, y) ** 2

    def T(x, y):
        K = 1
        temp = 1.0 / density(x, y) * (K - b2(x, y) * 0.5)
        assert np.all(temp > 0)
        return temp

    def vxyz(x, y):
        return 0.0

    def vthxyz(x, y):
        return np.sqrt(T(x, y))

    vHigh = {
        "vbulkx": vxyz,
        "vbulky": vxyz,
        "vbulkz": vxyz,
        "vthx": vthxyz,
        "vthy": vthxyz,
        "vthz": vthxyz,
        "nbr_part_per_cell": 90,
    }

    vLow = {
        "vbulkx": vxyz,
        "vbulky": vxyz,
        "vbulkz": vxyz,
        "vthx": vthxyz,
        "vthy": vthxyz,
        "vthz": vthxyz,
        "nbr_part_per_cell": 120,
    }

    ph.MaxwellianFluidModel(
        bx=bx,
        by=by,
        bz=bz,
        brotons={"charge": 1, "density": densityHigh, **vHigh, "init": {"seed": 12334}},
        protons={"charge": 1, "density": densityLow, **vLow, "init": {"seed": 43210}},
    )
    ph.ElectronModel(closure="isothermal", Te=0.0)

    for pop in ["brotons", "protons"]:
        ph.ParticleDiagnostics(
            quantity="domain", write_timestamps=timestamps, population_name=pop
        )

    if loadbalancing:
        ph.LoadBalancer(**loadbalancing)

    return sim


def get_time(path, time=0, pop="protons", datahier=None):
    time = "{:.10f}".format(time)
    from pyphare.pharesee.hierarchy import hierarchy_from

    return hierarchy_from(
        h5_filename=path + f"/ions_pop_{pop}_domain.h5", times=time, hier=datahier
    )


def get_particles(diag_dir, time=0):
    hier = get_time(diag_dir, time)
    hier = get_time(diag_dir, time, "brotons", hier)
    return hier


def time_info(diag_dir, time=0):
    hier = get_particles(diag_dir, time)

    per_rank = {f"p{rank}": 0 for rank in range(mpi_size)}

    def _parse_rank(patch_id):
        return patch_id.split("#")[0]

    for ilvl, lvl in hier.levels().items():
        for patch in lvl:
            for pd_key, pd in patch.patch_datas.items():
                per_rank[_parse_rank(patch.id)] += pd.size()

    return per_rank


@ddt
class LoadBalancingTest(SimulatorTest):
    def tearDown(self):
        ph.global_vars.sim = None

    def run_sim(self, diags_dir, dic={}):
        ph.global_vars.sim = None
        self.register_diag_dir_for_cleanup(diags_dir)
        Simulator(config(diags_dir, dic)).run()
        return diags_dir

    @data(dict(auto=True, every=1))
    @unpack
    def test_raises(self, **lbkwargs):
        if mpi_size == 1:  # doesn't make sense
            return

        with self.assertRaises(RuntimeError):
            diag_dir = self.run_sim(
                self.unique_diag_dir_for_test_case(diag_outputs, ndim, interp),
                dict(active=True, mode="nppc", tol=0.01, **lbkwargs),
            )
            # does not get here

    @unittest.skip("should change with moments")
    @data(
        dict(auto=True),  # tolerance checks
        dict(on_init=True, every=0),  # on init only
        dict(on_init=True, every=1),
        dict(on_init=False, auto=True, next_rebalance=1),
        dict(on_init=False, every=1),
    )
    @unpack
    def test_has_balanced(self, **lbkwargs):
        if mpi_size == 1:  # doesn't make sense
            return

        diag_dir = self.run_sim(
            self.unique_diag_dir_for_test_case(diag_outputs, ndim, interp),
            dict(active=True, mode="nppc", tol=0.01, **lbkwargs),
        )

        if cpp.mpi_rank() == 0:
            t0_sdev = np.std(list(time_info(diag_dir).values()))
            tend_sdev = np.std(list(time_info(diag_dir, timestamps[-1]).values()))
            self.assertLess(tend_sdev, t0_sdev * 0.15)  # empirical

    @unittest.skip("should change with moments")
    def test_has_not_balanced_as_defaults(self):
        if mpi_size == 1:  # doesn't make sense
            return

        diag_dir = self.run_sim(
            self.unique_diag_dir_for_test_case(diag_outputs, ndim, interp)
        )

        if cpp.mpi_rank() == 0:
            t0_sdev = np.std(list(time_info(diag_dir).values()))
            tend_sdev = np.std(list(time_info(diag_dir, timestamps[-1]).values()))
            self.assertGreater(tend_sdev, t0_sdev * 0.1)  # empirical

    @unittest.skip("should change with moments")
    def test_compare_is_and_is_not_balanced(self):
        if mpi_size == 1:  # doesn't make sense
            return

        check_time = 0.001

        not_hier = get_particles(
            self.run_sim(
                self.unique_diag_dir_for_test_case(diag_outputs + "/not", ndim, interp)
            ),
            check_time,
        )
        is_hier = get_particles(
            self.run_sim(
                self.unique_diag_dir_for_test_case(diag_outputs, ndim, interp),
                dict(active=True, auto=True, mode="nppc", tol=0.05),
            ),
            check_time,
        )

        if cpp.mpi_rank() == 0:
            not_hier = single_patch_per_level_per_pop_from(not_hier)
            is_hier = single_patch_per_level_per_pop_from(is_hier)
            self.assertTrue(hierarchy_compare(not_hier, is_hier))


if __name__ == "__main__":
    unittest.main()
