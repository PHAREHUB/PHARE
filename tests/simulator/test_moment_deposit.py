#!/usr/bin/env python3
"""
After a few time steps, the density and flux of each population must equal the deposit of
the domain particles of all patches of the level, whatever the order in which the solver
deposits, exchanges and sums them. The deposit is recomputed here with B-splines and compared
on every patch, on the nodes that level-ghost particles do not reach. Diagnostics are stored
in single precision, hence the tolerance; a missing or double-counted particle is ~1e-2.
"""

import unittest
import itertools
import numpy as np
from types import SimpleNamespace
from ddt import ddt, data, unpack

import pyphare.pharein as ph
from pyphare.pharesee.run import Run
from pyphare.simulator.simulator import Simulator

from tests.simulator import SimulatorTest

ph.NO_GUI()

# the beam drifts by half a cell over the run, so about 5% of the particles of a 10 cell wide
# patch cross one of its borders
time_step = 0.005
time_step_nbr = 10
final_time = time_step * time_step_nbr
pop_names = ["protons", "beam"]
cells = {1: (60,), 2: (30, 30), 3: (12, 12, 12)}
dl = {1: 0.2, 2: 0.2, 3: 0.3}

# nodes this far from a level border receive no level-ghost contribution:
# particle ghost width (hybrid_options.hpp) + stencil half-width + 1
margin = {1: 1 + 1 + 1, 2: 2 + 2 + 1, 3: 2 + 2 + 1}


def bspline(t, interp):
    """centered cardinal B-spline of degree `interp` (PHARE interpolation order)"""
    t = np.abs(t)
    if interp == 1:
        return np.where(t < 1, 1 - t, 0.0)
    if interp == 2:
        return np.where(
            t < 0.5, 0.75 - t**2, np.where(t < 1.5, 0.5 * (1.5 - t) ** 2, 0.0)
        )
    return np.where(
        t < 1, 2 / 3 - t**2 + t**3 / 2, np.where(t < 2, (2 - t) ** 3 / 6, 0.0)
    )


def deposit(particles, lower, shape, interp, periodic):
    """density and flux on the primal nodes lower..lower+shape-1 (AMR indexes)"""
    ndim = len(shape)
    density = np.zeros(shape)
    flux = [np.zeros(shape) for _ in range(3)]
    position = particles.iCells + particles.deltas
    base = np.floor(position).astype(int)

    for offset in itertools.product(range(-2, 3), repeat=ndim):
        node = base + np.asarray(offset)
        weight = particles.weights.reshape(-1).copy()
        for d in range(ndim):
            weight *= bspline(position[:, d] - node[:, d], interp)
        local = node - np.asarray(lower)
        if periodic:  # the last node of each direction is the image of the first
            local %= np.asarray(shape) - 1
        inside = np.all((local >= 0) & (local < np.asarray(shape)), axis=1) & (
            weight != 0
        )
        index = tuple(local[inside].T)
        np.add.at(density, index, weight[inside])
        for c in range(3):
            np.add.at(flux[c], index, weight[inside] * particles.v[inside, c])

    if periodic:
        for array in [density] + flux:
            for d in range(ndim):
                np.moveaxis(array, d, 0)[-1] = np.moveaxis(array, d, 0)[0]
    return density, flux


def config(ndim, interp, diag_dir):
    boxes = {
        1: {"L0": {"B0": [(20,), (40,)]}},
        2: {"L0": {"B0": [(8, 8), (21, 21)]}},
        3: {"L0": {"B0": [(3, 3, 3), (8, 8, 8)]}},
    }[ndim]
    sim = ph.Simulation(
        time_step=time_step,
        time_step_nbr=time_step_nbr,
        cells=cells[ndim],
        dl=[dl[ndim]] * ndim,
        interp_order=interp,
        refinement_boxes=boxes,
        largest_patch_size=8 if ndim == 3 else 12,
        smallest_patch_size=6,
        hyper_resistivity=0.002,
        strict=True,
        diag_options={
            "format": "phareh5",
            "options": {"dir": diag_dir, "mode": "overwrite"},
        },
    )
    L = [c * dl[ndim] for c in cells[ndim]]

    def zero(*xyz):
        return 0.0 * xyz[0]

    def one(*xyz):
        return 1.0 + 0.0 * xyz[0]

    def density(*xyz):
        return 1.0 + 0.3 * np.prod(
            [np.cos(2 * np.pi * x / l) for x, l in zip(xyz, L)], axis=0
        )

    def vx(*xyz):
        return 0.5 * np.sin(2 * np.pi * xyz[0] / L[0])

    def drift(*xyz):
        return 2.0 + 0.0 * xyz[0]

    def vth(*xyz):
        return 0.3 + 0.0 * xyz[0]

    def pop(dens, vbulkx, seed):
        return {
            "charge": 1,
            "density": dens,
            "vbulkx": vbulkx,
            "vbulky": zero,
            "vbulkz": zero,
            "vthx": vth,
            "vthy": vth,
            "vthz": vth,
            "nbr_part_per_cell": 10 if ndim < 3 else 4,
            "init": {"seed": seed},
        }

    ph.MaxwellianFluidModel(
        bx=one, by=zero, bz=zero, protons=pop(density, vx, 1), beam=pop(one, drift, 2)
    )
    ph.ElectronModel(closure="isothermal", Te=0.1)

    timestamps = [final_time]
    for name in pop_names:
        for quantity in ["density", "flux"]:
            ph.FluidDiagnostics(
                quantity=quantity, write_timestamps=timestamps, population_name=name
            )
        ph.ParticleDiagnostics(
            quantity="domain", write_timestamps=timestamps, population_name=name
        )
    return sim


@ddt
class MomentDepositTest(SimulatorTest):
    @data(*itertools.product([1, 2, 3], [1, 2, 3]))
    @unpack
    def test_moments_are_the_deposit_of_owned_particles(self, ndim, interp):
        diag_dir = f"phare_outputs/moment_deposit/{ndim}d_{interp}"
        self.register_diag_dir_for_cleanup(diag_dir)
        Simulator(config(ndim, interp, diag_dir)).run().reset()
        ph.global_vars.sim = None

        run = Run(diag_dir)
        crossed = 0

        for name in pop_names:
            particles = run.GetParticles(final_time, name)
            density = run.GetN(final_time, name)
            flux = run.GetFlux(final_time, name)

            crossed += self.crossings(
                particles.level(0, final_time).patches, name, ndim
            )

            for ilvl, level in particles.levels(final_time).items():
                datasets = [
                    p.patch_datas[f"{name}_domain"].dataset for p in level.patches
                ]
                merged = SimpleNamespace(
                    **{
                        attr: np.concatenate([getattr(ds, attr) for ds in datasets])
                        for attr in ["iCells", "deltas", "v", "weights"]
                    }
                )

                periodic = ilvl == 0
                lower = np.min([p.box.lower for p in level.patches], axis=0)
                upper = np.max([p.box.upper for p in level.patches], axis=0)
                if periodic:
                    lower = np.zeros(ndim, dtype=int)
                    upper = np.asarray(cells[ndim]) - 1
                shape = tuple(upper - lower + 2)
                expected_n, expected_f = deposit(merged, lower, shape, interp, periodic)

                for patch_n, patch_f in zip(
                    density.level(ilvl, final_time).patches,
                    flux.level(ilvl, final_time).patches,
                ):
                    box = patch_n.box
                    sl = tuple(
                        slice(box.lower[d] - lower[d], box.upper[d] - lower[d] + 2)
                        for d in range(ndim)
                    )
                    keep = np.ones([box.shape[d] + 1 for d in range(ndim)], dtype=bool)
                    if not periodic:
                        for d in range(ndim):
                            nodes = np.arange(box.lower[d], box.upper[d] + 2)
                            far = (nodes - lower[d] >= margin[interp]) & (
                                upper[d] + 1 - nodes >= margin[interp]
                            )
                            keep &= far.reshape(
                                [-1 if i == d else 1 for i in range(ndim)]
                            )
                    if not keep.any():
                        continue

                    actual = patch_n.patch_datas["value"][box]
                    np.testing.assert_allclose(
                        actual[keep], expected_n[sl][keep], rtol=1e-6, atol=1e-6
                    )
                    for c, comp in enumerate("xyz"):
                        actual = patch_f.patch_datas[comp][box]
                        np.testing.assert_allclose(
                            actual[keep], expected_f[c][sl][keep], rtol=1e-6, atol=1e-6
                        )

        print(f"{ndim}D interp {interp}: {crossed} particles crossed a patch border")
        # the test is meaningful only if particles changed patch during the run
        self.assertGreater(crossed, 0)

    @staticmethod
    def crossings(patches, name, ndim):
        """particles whose straight-line position at t=0 lies in another patch, along x"""
        crossed = 0
        for patch in patches:
            ds = patch.patch_datas[f"{name}_domain"].dataset
            back = (
                ds.iCells[:, 0] + ds.deltas[:, 0] - ds.v[:, 0] * final_time / dl[ndim]
            )
            crossed += np.count_nonzero(
                (back < patch.box.lower[0]) | (back > patch.box.upper[0] + 1)
            )
        return crossed


if __name__ == "__main__":
    unittest.main()
