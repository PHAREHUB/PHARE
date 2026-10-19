#!/usr/bin/env python3
"""
Ion-ion collisions in a full hybrid run: a uniform plasma along B with a temperature
anisotropy below the firehose and mirror thresholds relaxes at the NRL rate, populations of
one species collide as one, and a restart gives bitwise the run that did not stop.
"""

import copy
import unittest
import numpy as np
from ddt import ddt, data

import pyphare.pharein as ph
from pyphare import cpp
from pyphare.pharesee.run import Run
from pyphare.simulator.simulator import Simulator

from tests.simulator import SimulatorTest

ph.NO_GUI()

nu0 = 0.1
time_step = 0.002  # s at the thermal speed ~0.014, rates within ~3%
time_step_nbr = 110  # about a quarter of the anisotropy relaxes
final_time = time_step * time_step_nbr
vth_par, vth_perp = 0.2, 0.3  # beta_par 0.08, beta_perp 0.18
# 50k particles: the measured decay fluctuates by about 4.5% from one sampling to another
cells = 24
ppc = 2100


def config(
    diag_dir, collisions=True, populations=1, restart_options=None, refined=False
):
    # the refined level is made of split particles, of weights unequal within a fine cell
    sim = ph.Simulation(
        refinement_boxes={"L0": {"B0": [(6,), (17,)]}} if refined else None,
        time_step=time_step,
        time_step_nbr=time_step_nbr,
        cells=cells,
        dl=0.2,
        interp_order=1,
        largest_patch_size=6,
        smallest_patch_size=6,
        hyper_resistivity=0.001,
        collisions={"nu0": nu0, "seed": 1} if collisions else None,
        diag_options={
            "format": "phareh5",
            "options": {"dir": diag_dir, "mode": "overwrite"},
        },
        restart_options=restart_options,
    )

    def one(x):
        return 1.0 + 0.0 * x

    def zero(x):
        return 0.0 * x

    def pop(seed):
        return {
            "charge": 1,
            "mass": 1,
            "density": lambda x: one(x) / populations,
            "vbulkx": zero,
            "vbulky": zero,
            "vbulkz": zero,
            "vthx": lambda x: vth_par + zero(x),
            "vthy": lambda x: vth_perp + zero(x),
            "vthz": lambda x: vth_perp + zero(x),
            "nbr_part_per_cell": ppc // populations,
            "init": {"seed": seed},
        }

    ph.MaxwellianFluidModel(
        bx=one,
        by=zero,
        bz=zero,
        **{f"protons{i}": pop(i + 1) for i in range(populations)},
    )
    ph.ElectronModel(closure="isothermal", Te=0.04)
    restarted = restart_options is not None and "restart_time" in restart_options
    timestamps = [final_time] if restarted else [0.0, final_time]
    for name in [f"protons{i}" for i in range(populations)]:
        ph.ParticleDiagnostics(
            quantity="domain", write_timestamps=timestamps, population_name=name
        )
    return sim


def temperatures(run, time, populations, level=0):
    """T_par, T_perp and density of all the domain particles of a level"""
    v, w, level_cells = [], [], 0
    for i in range(populations):
        for patch in run.GetParticles(time, f"protons{i}").level(level, time).patches:
            particles = patch.patch_datas[f"protons{i}_domain"].dataset
            v.append(particles.v)
            w.append(particles.weights.reshape(-1))
            level_cells += patch.box.shape[0] if i == 0 else 0
    v, w = np.concatenate(v), np.concatenate(w)
    mean = np.average(v, axis=0, weights=w)
    T = np.average((v - mean) ** 2, axis=0, weights=w)
    return T[0], 0.5 * (T[1] + T[2]), w.sum() / level_cells


def nrl_decay(T_par, T_perp, density, t):
    """ln of the anisotropy left after t, integrating the NRL isotropization rate"""

    def rate(T):
        A = T[1] / T[0] - 1.0
        r = np.sqrt(abs(A))
        g = np.arctan(r) / r if A > 0 else np.arctanh(r) / r
        nu_T = (
            nu0
            * density
            / (2 * np.sqrt(np.pi) * T[0] ** 1.5)
            * (-3 + (A + 3) * g)
            / A**2
        )
        return np.array([2 * nu_T * (T[1] - T[0]), -nu_T * (T[1] - T[0])])

    T, steps = np.array([T_par, T_perp]), 2000
    for _ in range(steps):
        k1 = rate(T)
        k2 = rate(T + 0.5 * t / steps * k1)
        k3 = rate(T + 0.5 * t / steps * k2)
        k4 = rate(T + t / steps * k3)
        T = T + t / steps / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    return np.log((T[1] - T[0]) / (T_perp - T_par))


@ddt
class CollisionsTest(SimulatorTest):
    def run_sim(self, diag_dir, **kwargs):
        self.register_diag_dir_for_cleanup(diag_dir)
        Simulator(config(diag_dir, **kwargs)).run().reset()
        ph.global_vars.sim = None
        return Run(diag_dir)

    def decay_ratio(self, run, populations, level=0):
        T_par0, T_perp0, density = temperatures(run, 0.0, populations, level)
        T_par, T_perp, _ = temperatures(run, final_time, populations, level)
        simulated = np.log((T_perp - T_par) / (T_perp0 - T_par0))
        return simulated / nrl_decay(T_par0, T_perp0, density, final_time)

    @data(1, 2)
    def test_anisotropy_relaxes_at_the_nrl_rate(self, populations):
        # two populations of one species collide as one: same rate as a single population
        run = self.run_sim(
            f"phare_outputs/collisions/relax_{populations}", populations=populations
        )
        if cpp.mpi_rank() == 0:
            ratio = self.decay_ratio(run, populations)
            print(f"{populations} population(s): anisotropy decay over NRL {ratio:.3f}")
            self.assertAlmostEqual(ratio, 1.0, delta=0.15)

    def test_refined_level_relaxes_at_the_nrl_rate(self):
        run = self.run_sim("phare_outputs/collisions/refined", refined=True)
        if cpp.mpi_rank() == 0:
            ratio = self.decay_ratio(run, 1, level=1)
            print(f"refined level: anisotropy decay over NRL {ratio:.3f}")
            self.assertAlmostEqual(ratio, 1.0, delta=0.15)

    def test_refined_level_without_collisions_keeps_its_anisotropy(self):
        # split particles of unequal weights must not relax by themselves (Lehe et al. 2026
        # warn that bare PIC drives different weights towards equal weight x temperature)
        run = self.run_sim(
            "phare_outputs/collisions/refined_off", collisions=False, refined=True
        )
        if cpp.mpi_rank() == 0:
            T_par0, T_perp0, _ = temperatures(run, 0.0, 1, level=1)
            T_par, T_perp, _ = temperatures(run, final_time, 1, level=1)
            self.assertAlmostEqual(
                (T_perp - T_par) / (T_perp0 - T_par0), 1.0, delta=0.03
            )

    def test_no_collisions_no_relaxation(self):
        run = self.run_sim("phare_outputs/collisions/off", collisions=False)
        if cpp.mpi_rank() == 0:
            T_par0, T_perp0, _ = temperatures(run, 0.0, 1)
            T_par, T_perp, _ = temperatures(run, final_time, 1)
            self.assertAlmostEqual(
                (T_perp - T_par) / (T_perp0 - T_par0), 1.0, delta=0.03
            )

    def test_restart_gives_the_same_particles(self):
        if cpp.mpi_size() > 1:
            # restarts are not bitwise under MPI without collisions either, already on master
            # (ba151d929, 3 ranks: 2221 of 50400 particles differ); cause not identified
            self.skipTest("bitwise restarts are serial only")
        restart_time = time_step * (time_step_nbr // 2)
        full_dir = "phare_outputs/collisions/restart_full"
        restarts = {"dir": full_dir, "mode": "overwrite", "timestamps": [restart_time]}
        full = self.run_sim(full_dir, restart_options=copy.deepcopy(restarts))

        restarted_dir = "phare_outputs/collisions/restart_restarted"
        restarts["restart_time"] = restart_time
        del restarts["timestamps"]
        restarted = self.run_sim(restarted_dir, restart_options=restarts)

        if cpp.mpi_rank() > 0:
            return
        for run0, run1 in [(full, restarted)]:
            patches0 = (
                run0.GetParticles(final_time, "protons0").level(0, final_time).patches
            )
            patches1 = (
                run1.GetParticles(final_time, "protons0").level(0, final_time).patches
            )
            for p0, p1 in zip(patches0, patches1):
                self.assertEqual(p0.box, p1.box)
                d0, d1 = (p.patch_datas["protons0_domain"].dataset for p in (p0, p1))

                def ordered(d):
                    order = np.lexsort((d.deltas[:, 0], d.iCells[:, 0]))
                    return d.v[order], d.weights[order]

                v0, w0 = ordered(d0)
                v1, w1 = ordered(d1)
                np.testing.assert_array_equal(v0, v1)
                np.testing.assert_array_equal(w0, w1)


if __name__ == "__main__":
    unittest.main()
