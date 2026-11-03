import itertools
import unittest
from copy import deepcopy

import numpy as np
from ddt import data, ddt, unpack

import pyphare.pharein as ph
from pyphare import cpp
from pyphare.core import phare_utilities as phut
from pyphare.simulator.simulator import Simulator, startMPI
from tests.diagnostic import all_timestamps, dump_all_diags
from tests.simulator import SimulatorTest

# only in 2d for now
ppc_per_dim = [100, 25, 10]


def config(sim):
    ppc = ppc_per_dim[sim.ndim - 1]

    L = 0.5

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
    }

    ph.MaxwellianFluidModel(
        bx=bx,
        by=by,
        bz=bz,
        protons={
            "mass": 1,
            "charge": 1,
            "density": density,
            **vvv,
            "nbr_part_per_cell": ppc,
            "init": {"seed": 1337},
        },
        alpha={
            "mass": 4,
            "charge": 2,
            "density": density,
            **vvv,
            "nbr_part_per_cell": ppc,
            "init": {"seed": 2334},
        },
    )
    ph.ElectronModel(closure="isothermal", Te=0.12)
    return sim


simArgs = {
    "time_step_nbr": 1,
    "final_time": 0.001,
    "cells": 40,
    "dl": 0.3,
    "diag_options": {
        "format": "pharevtkhdf",
        "options": {"mode": "overwrite", "dir": "phare_outputs/vtk_diagnostic_test"},
    },
}


TENSOR_COMPONENTS = ("xx", "xy", "xz", "yy", "yz", "zz")


def vtkhdf_patches(path, time):
    """
    {ilvl: [(amr_box_lower, {name: primal nodes})]} at `time`, read with h5py only.
    Rows are per box in AMRBox order, nodes x fastest, each box repeated
    2 (2d) or 4 (1d) times to fill a 3d cell; the column count says scalar,
    vector (x,y,z) or symmetric tensor (xx,xy,xz,yy,yz,zz).
    """
    import h5py

    names = {1: ("",), 3: ("x", "y", "z"), 6: TENSOR_COMPONENTS}
    patches = {}
    with h5py.File(path, "r") as f:
        root = f["VTKHDF"]
        step = int(np.argmin(np.abs(root["Steps/Values"][:] - time)))
        assert np.isclose(root["Steps/Values"][step], time)

        ilvl = 0
        while f"Level{ilvl}" in root:
            steps = root[f"Steps/Level{ilvl}"]
            box0 = int(steps["AMRBoxOffset"][step])
            nbox = int(steps["NumberOfAMRBox"][step])
            row = int(steps["PointDataOffset/data"][step])
            data = root[f"Level{ilvl}/PointData/data"]
            boxes = root[f"Level{ilvl}/AMRBox"][box0 : box0 + nbox]
            ndim = len(root[f"Level{ilvl}"].attrs["Spacing"])
            x_times = [4, 2, 1][ndim - 1]

            patches[ilvl] = []
            for box in boxes:
                lower = box[0 : 2 * ndim : 2]
                shape = box[1 : 2 * ndim : 2] - lower + 2  # primal nodes
                size = int(np.prod(shape))
                rows = data[row : row + size]
                row += size * x_times
                patches[ilvl].append(
                    (
                        tuple(lower),
                        {
                            name: rows[:, c].reshape(shape, order="F")
                            for c, name in enumerate(names[data.shape[1]])
                        },
                    )
                )
            ilvl += 1
    return patches


def phareh5_patches(path, time, quantity):
    """same as vtkhdf_patches for a phareh5 file, ghost nodes dropped"""
    import h5py

    patches = {}
    with h5py.File(path, "r") as f:
        times = {float(t): t for t in f["t"]}
        time_grp = f["t"][times[min(times, key=lambda t: abs(t - time))]]
        for lvl_key, lvl in time_grp.items():
            ilvl = int(lvl_key[2:])  # pl1 -> 1
            patches[ilvl] = []
            for patch in lvl.values():
                lower, upper = patch.attrs["lower"], patch.attrs["upper"]
                shape = upper - lower + 2  # primal nodes
                fields = {}
                for name in TENSOR_COMPONENTS:
                    ds = patch[f"{quantity}_{name}"]
                    g = int(ds.attrs["ghosts"])
                    fields[name] = ds[tuple(slice(g, g + n) for n in shape)]
                patches[ilvl].append((tuple(lower), fields))
    return patches


def uniform_config(sim, density, vbulk, vth):
    """a uniform drifting Maxwellian, its pressure tensor is n m diag(vth**2)"""

    def const(value):
        return lambda *xyz: np.full(np.asarray(xyz[0]).shape, value)

    ph.MaxwellianFluidModel(
        bx=const(1.0),
        by=const(0.0),
        bz=const(0.0),
        protons={
            "mass": 1,
            "charge": 1,
            "density": const(density),
            **{f"vbulk{d}": const(v) for d, v in zip("xyz", vbulk)},
            **{f"vth{d}": const(v) for d, v in zip("xyz", vth)},
            "nbr_part_per_cell": 100,
            "init": {"seed": 1337},
        },
    )
    ph.ElectronModel(closure="isothermal", Te=0.12)
    return sim


def permute(dic):
    ndims = [2]
    interp_orders = [1]
    dic.update(simArgs.copy())
    return [
        dict(
            ndim=ndim,
            interp=interp_order,
            simInput=deepcopy(dic),
        )
        for ndim, interp_order in itertools.product(ndims, interp_orders)
    ]


@ddt
class VTKDiagnosticsTest(SimulatorTest):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.simulator = None

    def _run(self, ndim, interp, simInput, diag_dir="", **kwargs):
        for key in ["cells", "dl"]:
            simInput[key] = list(phut.np_array_ify(simInput[key], ndim))
        simulation = config(self.simulation(interp_order=interp, **simInput))
        self.assertTrue(len(simulation.cells) == ndim)
        dump_all_diags(simulation.model.populations)
        Simulator(simulation).run().reset()
        return simulation.diag_options["options"]["dir"]

    @data(*permute({}))
    @unpack
    def test_dump_diags(self, ndim, interp, simInput):
        print(f"test_dump_diags dim/interp:{ndim}/{interp}")

        b0 = [[10 for i in range(ndim)], [19 for i in range(ndim)]]
        simInput["refinement_boxes"] = {"L0": {"B0": b0}}
        local_out = self._run(ndim, interp, simInput)

        if cpp.mpi_rank() == 0:
            import h5py

            for pop, mass, charge in [("protons", 1, 1), ("alpha", 4, 2)]:
                for qty in ["density", "charge_density", "flux"]:
                    file = f"{local_out}/ions_pop_{pop}_{qty}.vtkhdf"
                    with h5py.File(file, "r") as h5_file:
                        self.assertEqual(h5_file.attrs["pop_mass"], mass)
                        self.assertEqual(h5_file.attrs["pop_charge"], charge)

            momentum_tensor_files = [f"{local_out}/ions_momentum_tensor.vtkhdf"] + [
                f"{local_out}/ions_pop_{pop}_momentum_tensor.vtkhdf"
                for pop in ["protons", "alpha"]
            ]
            for file in momentum_tensor_files:
                with h5py.File(file, "r") as h5_file:
                    data = h5_file["VTKHDF/Level0/PointData/data"]
                    self.assertEqual(data.shape[1], len(TENSOR_COMPONENTS))
                    self.assertGreater(data.shape[0], 0)

            try:
                from pyphare.pharesee.phare_vtk import plot as plot_vtk

                plot_vtk(local_out + "/EM_B.vtkhdf", f"B{ndim}d_interp{interp}.vtk.png")
                plot_vtk(local_out + "/EM_E.vtkhdf", f"E{ndim}d_interp{interp}.vtk.png")
            except ModuleNotFoundError:
                print("WARNING: vtk python module not found - cannot make plots")

    def _run_momentum_tensor(self, ndim, interp, simInput, fmt):
        simInput = deepcopy(simInput)
        for key in ["cells", "dl"]:
            simInput[key] = list(phut.np_array_ify(simInput[key], ndim))
        b0 = [[10 for i in range(ndim)], [19 for i in range(ndim)]]
        simInput["refinement_boxes"] = {"L0": {"B0": b0}}
        simInput["diag_options"]["format"] = fmt
        simInput["diag_options"]["options"]["dir"] += f"/{fmt}"
        simulation = config(self.simulation(interp_order=interp, **simInput))
        timestamps = all_timestamps(simulation)
        ph.FluidDiagnostics(quantity="momentum_tensor", write_timestamps=timestamps)
        for pop in simulation.model.populations:
            ph.FluidDiagnostics(
                quantity="momentum_tensor",
                write_timestamps=timestamps,
                population_name=pop,
            )
        Simulator(simulation).run().reset()
        return simulation.diag_options["options"]["dir"], timestamps

    @data(*permute({}))
    @unpack
    def test_momentum_tensor_matches_phareh5(self, ndim, interp, simInput):
        """same run dumped in both formats: the tensors must match on every level"""
        h5_out, timestamps = self._run_momentum_tensor(
            ndim, interp, simInput, "phareh5"
        )
        vtk_out, _ = self._run_momentum_tensor(ndim, interp, simInput, "pharevtkhdf")

        if cpp.mpi_rank() == 0:
            for basename in [
                "ions_momentum_tensor",
                "ions_pop_protons_momentum_tensor",
                "ions_pop_alpha_momentum_tensor",
            ]:
                for time in timestamps:
                    h5 = phareh5_patches(
                        f"{h5_out}/{basename}.h5", time, "momentum_tensor"
                    )
                    vtk = vtkhdf_patches(f"{vtk_out}/{basename}.vtkhdf", time)
                    self.assertEqual(sorted(h5), sorted(vtk))
                    self.assertGreater(len(vtk), 1)  # the refined level is there
                    for ilvl in vtk:
                        h5_lvl = dict(h5[ilvl])
                        self.assertEqual(
                            sorted(h5_lvl), sorted(l for l, _ in vtk[ilvl])
                        )
                        for lower, fields in vtk[ilvl]:
                            for name in TENSOR_COMPONENTS:
                                vtk_data = fields[name]
                                tol = max(1e-10, 8 * np.finfo(vtk_data.dtype).eps)
                                # under MPI two identical runs only agree to rounding:
                                #  particle and ghost-sum order follow message arrival
                                np.testing.assert_allclose(
                                    vtk_data,
                                    h5_lvl[lower][name],
                                    rtol=tol,
                                    atol=tol,
                                    err_msg=f"{basename} M{name} t={time} L{ilvl}",
                                )

    @data(*permute({}))
    @unpack
    def test_pressure_tensor_of_uniform_maxwellian(self, ndim, interp, simInput):
        """P = M - rho V V formed from the vtkhdf files, against n m diag(vth**2)"""
        density, vbulk, vth = 1.0, (0.2, -0.1, 0.3), (0.3, 0.4, 0.5)

        for key in ["cells", "dl"]:
            simInput[key] = list(phut.np_array_ify(simInput[key], ndim))
        b0 = [[10 for i in range(ndim)], [19 for i in range(ndim)]]
        simInput["refinement_boxes"] = {"L0": {"B0": b0}}
        simulation = uniform_config(
            self.simulation(interp_order=interp, **simInput), density, vbulk, vth
        )
        timestamps = all_timestamps(simulation)
        for quantity in ["mass_density", "bulkVelocity", "momentum_tensor"]:
            ph.FluidDiagnostics(quantity=quantity, write_timestamps=timestamps)
        Simulator(simulation).run().reset()
        local_out = simulation.diag_options["options"]["dir"]

        if cpp.mpi_rank() == 0:
            time = timestamps[-1]
            rho = vtkhdf_patches(f"{local_out}/ions_mass_density.vtkhdf", time)
            V = vtkhdf_patches(f"{local_out}/ions_bulkVelocity.vtkhdf", time)
            M = vtkhdf_patches(f"{local_out}/ions_momentum_tensor.vtkhdf", time)

            expected = density * np.diag(np.asarray(vth) ** 2)
            index = {"x": 0, "y": 1, "z": 2}
            for ilvl in M:
                P = {name: [] for name in TENSOR_COMPONENTS}
                for (_, r), (_, v), (_, m) in zip(rho[ilvl], V[ilvl], M[ilvl]):
                    for name in TENSOR_COMPONENTS:
                        i, j = name
                        P[name] += [(m[name] - r[""] * v[i] * v[j]).ravel()]

                for name in TENSOR_COMPONENTS:
                    i, j = index[name[0]], index[name[1]]
                    mean = np.mean(np.concatenate(P[name]))
                    # >= 4e4 particles per level: sampling noise on the mean is below 1%
                    self.assertAlmostEqual(
                        mean,
                        expected[i, j],
                        delta=0.03 * expected.max(),
                        msg=f"P{name} L{ilvl}",
                    )

    @data(*permute({}))
    @unpack
    def test_unsupported_quantity_raises(self, ndim, interp, simInput):
        """a quantity the vtkhdf fluid writer cannot write must not be skipped silently"""
        for key in ["cells", "dl"]:
            simInput[key] = list(phut.np_array_ify(simInput[key], ndim))
        simulation = config(self.simulation(interp_order=interp, **simInput))
        ph.FluidDiagnostics(  # per population mass density: no vtkhdf writer for it
            quantity="mass_density",
            write_timestamps=all_timestamps(simulation),
            population_name="protons",
        )
        simulator = Simulator(simulation)
        try:
            with self.assertRaises(ValueError) as error:
                simulator.initialize()
            self.assertIn(
                "unsupported fluid diagnostic", str(error.exception.__cause__)
            )
        finally:
            simulator.reset()


if __name__ == "__main__":
    startMPI()
    unittest.main()
