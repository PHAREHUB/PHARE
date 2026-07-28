#!/usr/bin/env python3

import numpy as np
import pyphare.pharein as ph
from pyphare.pharesee.hierarchy.fromh5 import get_all_available_quantities_from_h5
from pyphare.pharesee.hierarchy.patchdata import ParticleData
from pyphare.simulator.simulator import Simulator, startMPI

from pyphare import cpp
from tests.simulator import SimulatorTest

ph.NO_GUI()


time_step = 0.005
time_step_nbr = 4
final_time = time_step * time_step_nbr

# ModelView::MTAlgo::getOrCreateSchedule (src/diagnostic/diagnostic_model_view.hpp)
# caches the border-overlap schedule it uses to compute the (per population) momentum/
# pressure tensor by SAMRAI level number only, and never rebuilds it when that level is
# regridded (unlike the regular messenger pipeline, see multiphysics_integrator.hpp).
# L2 is regridded by tagging on essentially every substep of its parent level, so a schedule
# built at t=0 is expected to already be stale relative to L2's actual box layout well before
# `final_time`. Two candidate restart points are used below so this does not depend on
# precisely predicting SAMRAI's regrid cadence.
restart_times = [time_step, 3 * time_step]  # 0.005, 0.015

diag_dir = "phare_outputs/test_regridding"


def config(diag_dir_, write_timestamps, restart_source_dir=None, restart_time=None):
    L = 0.5

    restart_options = {"dir": restart_source_dir or diag_dir_, "mode": "overwrite"}
    if restart_source_dir is None:
        restart_options["timestamps"] = restart_times
    if restart_time is not None:
        restart_options["restart_time"] = restart_time

    sim = ph.Simulation(
        time_step=time_step,
        time_step_nbr=time_step_nbr,
        cells=(40, 40),
        dl=(0.40, 0.40),
        refinement="tagging",
        max_nbr_levels=3,
        nesting_buffer=1,
        clustering="tile",
        tag_buffer="1",
        hyper_resistivity=0.002,
        resistivity=0.001,
        diag_options={
            "format": "phareh5",
            "options": {"dir": diag_dir_, "mode": "overwrite"},
        },
        restart_options=restart_options,
    )

    def density(x, y):
        Ly = sim.simulation_domain()[1]
        return (
            0.4
            + 1.0 / np.cosh((y - Ly * 0.3) / L) ** 2
            + 1.0 / np.cosh((y - Ly * 0.7) / L) ** 2
        )

    def S(y, y0, l):
        return 0.5 * (1.0 + np.tanh((y - y0) / l))

    def by(x, y):
        Lx = sim.simulation_domain()[0]
        Ly = sim.simulation_domain()[1]
        sigma = 1.0
        dB = 0.1

        x0 = x - 0.5 * Lx
        y1 = y - 0.3 * Ly
        y2 = y - 0.7 * Ly

        dBy1 = 2 * dB * x0 * np.exp(-(x0**2 + y1**2) / (sigma) ** 2)
        dBy2 = -2 * dB * x0 * np.exp(-(x0**2 + y2**2) / (sigma) ** 2)

        return dBy1 + dBy2

    def bx(x, y):
        Lx = sim.simulation_domain()[0]
        Ly = sim.simulation_domain()[1]
        sigma = 1.0
        dB = 0.1

        x0 = x - 0.5 * Lx
        y1 = y - 0.3 * Ly
        y2 = y - 0.7 * Ly

        dBx1 = -2 * dB * y1 * np.exp(-(x0**2 + y1**2) / (sigma) ** 2)
        dBx2 = 2 * dB * y2 * np.exp(-(x0**2 + y2**2) / (sigma) ** 2)

        v1 = -1
        v2 = 1.0
        return v1 + (v2 - v1) * (S(y, Ly * 0.3, L) - S(y, Ly * 0.7, L)) + dBx1 + dBx2

    def bz(x, y):
        return 0.0

    def b2(x, y):
        return bx(x, y) ** 2 + by(x, y) ** 2 + bz(x, y) ** 2

    def T(x, y):
        K = 0.7
        temp = 1.0 / density(x, y) * (K - b2(x, y) * 0.5)
        assert np.all(temp > 0)
        return temp

    def vx(x, y):
        return 0.0

    def vy(x, y):
        return 0.0

    def vz(x, y):
        return 0.0

    def vthx(x, y):
        return np.sqrt(T(x, y))

    def vthy(x, y):
        return np.sqrt(T(x, y))

    def vthz(x, y):
        return np.sqrt(T(x, y))

    vvv = {
        "vbulkx": vx,
        "vbulky": vy,
        "vbulkz": vz,
        "vthx": vthx,
        "vthy": vthy,
        "vthz": vthz,
        "nbr_part_per_cell": 44,
    }

    ph.MaxwellianFluidModel(
        bx=bx,
        by=by,
        bz=bz,
        protons={
            "charge": 1,
            "density": density,
            **vvv,
        },
    )
    ph.ElectronModel(closure="isothermal", Te=0.0)

    for quantity in ["E", "B"]:
        ph.ElectromagDiagnostics(quantity=quantity, write_timestamps=write_timestamps)

    for quantity in [
        "mass_density",
        "charge_density",
        "bulkVelocity",
        "pressure_tensor",
    ]:
        ph.FluidDiagnostics(quantity=quantity, write_timestamps=write_timestamps)

    pop = "protons"
    ph.ParticleDiagnostics(
        quantity="domain", write_timestamps=write_timestamps, population_name=pop
    )

    for quantity in ["density", "charge_density", "pressure_tensor"]:
        ph.FluidDiagnostics(
            quantity=quantity, write_timestamps=write_timestamps, population_name=pop
        )
    return sim


class RegriddingTest(SimulatorTest):
    def __init__(self, *args, **kwargs):
        super(RegriddingTest, self).__init__(*args, **kwargs)
        self.simulator = None

    def tearDown(self):
        super(RegriddingTest, self).tearDown()
        if self.simulator is not None:
            self.simulator.reset()
        self.simulator = None
        ph.global_vars.sim = None

    def _run_continuous(self):
        diag_dir_a = f"{diag_dir}_continuous"
        run_timestamps = sorted({0.0, *restart_times, final_time})

        ph.global_vars.sim = None
        sim = config(diag_dir_a, run_timestamps)
        self.register_diag_dir_for_cleanup(diag_dir_a)
        Simulator(sim).run().reset()
        ph.global_vars.sim = None
        return diag_dir_a

    def _run_restarted(self, source_dir, restart_time):
        diag_dir_b = f"{diag_dir}_restart_{restart_time}"

        ph.global_vars.sim = None
        sim = config(
            diag_dir_b,
            [restart_time],
            restart_source_dir=source_dir,
            restart_time=restart_time,
        )
        self.register_diag_dir_for_cleanup(diag_dir_b)
        # only the state right at/after restore is needed: MTAlgo's schedule for this
        # process is built fresh on that first post-restart dump, so it cannot be stale.
        Simulator(sim).initialize().reset()
        ph.global_vars.sim = None
        return diag_dir_b

    def _assert_same_layout(self, hier0, hier1, time):
        self.assertEqual(
            len(hier0.levels()),
            len(hier1.levels()),
            f"level count mismatch at t={time}",
        )
        for ilvl in range(len(hier0.levels())):
            patches0 = hier0.level(ilvl).patches
            patches1 = hier1.level(ilvl).patches
            self.assertEqual(
                len(patches0),
                len(patches1),
                f"patch count mismatch at t={time}, level {ilvl}",
            )
            for patch0, patch1 in zip(patches0, patches1):
                self.assertEqual(
                    patch0.box, patch1.box, f"box mismatch at t={time}, level {ilvl}"
                )

    def _assert_same_data(self, hier0, hier1, time):
        for ilvl, lvl0 in hier0.levels().items():
            lvl1 = hier1.levels()[ilvl]
            for patch0, patch1 in zip(lvl0.patches, lvl1.patches):
                for key, pd0 in patch0.patch_datas.items():
                    pd1 = patch1.patch_datas[key]
                    if isinstance(pd0, ParticleData):
                        # particle ordering across identical physical state is not
                        # guaranteed to match, this isn't what this test is about
                        continue
                    np.testing.assert_equal(
                        pd0.dataset[:],
                        pd1.dataset[:],
                        err_msg=f"mismatch for '{key}' at t={time}, level {ilvl}",
                    )

    def test_regrid_does_not_go_stale(self):
        """
        Compares a continuous run against runs restarted right at candidate regrid
        points. Both reach the same particle/field state and the same patch layout at
        the comparison time, so any dataset mismatch between them is not physical - it
        exposes ModelView keeping a schedule built against a box layout that no longer
        exists.
        """
        diag_dir_a = self._run_continuous()
        restarted_dirs = [self._run_restarted(diag_dir_a, rt) for rt in restart_times]

        if cpp.mpi_rank() != 0:
            return

        for rt, diag_dir_b in zip(restart_times, restarted_dirs):
            with self.subTest(restart_time=rt):
                hier_a = get_all_available_quantities_from_h5(diag_dir_a, rt)
                hier_b = get_all_available_quantities_from_h5(diag_dir_b, rt)

                self.assertEqual(set(hier_a.quantities()), set(hier_b.quantities()))
                self._assert_same_layout(hier_a, hier_b, rt)
                self._assert_same_data(hier_a, hier_b, rt)


if __name__ == "__main__":
    import unittest

    startMPI()
    unittest.main()
