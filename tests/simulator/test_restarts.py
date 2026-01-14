#
#

import copy
import datetime
import unittest
import numpy as np
from time import sleep
from pathlib import Path

from datetime import timedelta
from ddt import ddt, data, unpack

import pyphare.pharein as ph

from pyphare.cpp import cpp_lib
from pyphare.pharesee.run import Run
from pyphare.simulator.simulator import Simulator

from tests.simulator import SimulatorTest
from tests.diagnostic import dump_all_diags
from pyphare.pharesee.hierarchy.patchdata import ParticleData
from pyphare.pharesee.hierarchy.fromh5 import get_all_available_quantities_from_h5

cpp = cpp_lib()


def permute(dic, expected_num_levels):
    # from pyphare.pharein.simulation import supported_dimensions # eventually
    dims = [1]  # supported_dimensions()
    return [
        [dim, interp, dic, expected_num_levels] for dim in dims for interp in [1, 2, 3]
    ]


def setup_model(ppc=100):
    def density(x):
        return 1.0

    def bx(x):
        return 0.0

    def S(x, x0, l):
        return 0.5 * (1 + np.tanh((x - x0) / l))

    def by(x):
        L = ph.global_vars.sim.simulation_domain()[0]
        v1, v2 = -1, 1.0
        return v1 + (v2 - v1) * (S(x, L * 0.25, 1) - S(x, L * 0.75, 1))

    def bz(x):
        return 0.5

    def b2(x):
        return bx(x) ** 2 + by(x) ** 2 + bz(x) ** 2

    def T(x):
        K = 1
        return 1 / density(x) * (K - b2(x) * 0.5)

    def vx(x):
        return 2.0

    def vy(x):
        return 0.0

    def vz(x):
        return 0.0

    def vxalpha(x):
        return 3.0

    def vthxyz(x):
        return T(x)

    vvv = {
        "vbulkx": vx,
        "vbulky": vy,
        "vbulkz": vz,
        "vthx": vthxyz,
        "vthy": vthxyz,
        "vthz": vthxyz,
    }
    vvvalpha = {
        "vbulkx": vxalpha,
        "vbulky": vy,
        "vbulkz": vz,
        "vthx": vthxyz,
        "vthy": vthxyz,
        "vthz": vthxyz,
    }
    model = ph.MaxwellianFluidModel(
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
            "mass": 4.0,
            "charge": 1,
            "density": density,
            **vvvalpha,
            "nbr_part_per_cell": ppc,
            "init": {"seed": 2334},
        },
    )
    ph.ElectronModel(closure="isothermal", Te=0.12)
    return model


timestep = 0.001
out = "phare_outputs/restarts"
simArgs = dict(
    # we are saving at timestep 4, and we have seen that restarted simulations with refinement boxes
    #  have regridding in places that don't exist in the original simulation
    #   we compare the immediate next timestep of both simulations with refinement boxes, as we have seen
    #   in this way neither simulations have any regrids, so are still comparable
    time_step_nbr=5,  # avoid regrid for refinement boxes https://github.com/LLNL/SAMRAI/issues/199
    time_step=timestep,
    boundary_types="periodic",
    cells=200,
    dl=0.3,
    diag_options=dict(format="phareh5", options=dict(dir=out, mode="overwrite")),
    restart_options=dict(dir=out, mode="overwrite"),
)


def dup(dic={}):
    args = copy.deepcopy(simArgs)
    args.update(dic)
    return args


@ddt
class RestartsTest(SimulatorTest):
    def __init__(self, *args, **kwargs):
        super(RestartsTest, self).__init__(*args, **kwargs)
        self.simulator = None

    def tearDown(self):
        super(RestartsTest, self).tearDown()
        if self.simulator is not None:
            self.simulator.reset()
        self.simulator = None
        ph.global_vars.sim = None

    def ddt_test_id(self):
        return self._testMethodName.split("_")[-1]

    def check_diags(self, diag_dir0, diag_dir1, pops, timestamps, expected_num_levels):
        if cpp.mpi_rank() > 0:
            return

        def count_levels_and_patches(qty):
            n_levels = len(qty.levels())
            n_patches = 0
            for ilvl in qty.levels().keys():
                n_patches += len(qty.level(ilvl).patches)
            return n_levels, n_patches

        self.assertGreater(len(timestamps), 0)
        for time in timestamps:
            checks = 0

            run0 = Run(diag_dir0)

            datahier0 = get_all_available_quantities_from_h5(diag_dir0, time)
            datahier1 = get_all_available_quantities_from_h5(diag_dir1, time)

            self.assertEqual(
                set(datahier0.quantities()),
                set(datahier1.quantities()),
            )

            self.assertEqual(len(datahier0.levels()), len(datahier1.levels()))
            for ilvl in range(len(datahier0.levels())):
                self.assertEqual(
                    len(datahier0.level(ilvl).patches),
                    len(datahier1.level(ilvl).patches),
                )
                for patch0, patch1 in zip(
                    datahier0.level(ilvl).patches, datahier1.level(ilvl).patches
                ):
                    self.assertEqual(patch0.box, patch1.box)

            self.assertGreater(len(datahier0.levels()), 0)

            for ilvl, lvl0 in datahier0.levels().items():
                patch_level1 = datahier1.levels()[ilvl]
                for p_idx, patch0 in enumerate(lvl0):
                    patch1 = patch_level1.patches[p_idx]
                    for pd_key, pd0 in patch0.patch_datas.items():
                        pd1 = patch1.patch_datas[pd_key]
                        self.assertNotEqual(id(pd0), id(pd1))
                        self.assertEqual(type(pd0), type(pd1))
                        if isinstance(pd1, ParticleData):
                            try:
                                self.assertEqual(pd0.dataset, pd1.dataset)
                            except AssertionError:
                                print(
                                    f"FAILED domain particles at time {time} {ilvl} {patch1.box} {patch0.box}"
                                )
                        else:
                            np.testing.assert_equal(pd0.dataset[:], pd1.dataset[:])
                        checks += 1

            n_levels, n_patches = count_levels_and_patches(
                run0.GetB(time, all_primal=False)
            )
            self.assertEqual(n_levels, expected_num_levels)
            self.assertGreaterEqual(n_patches, n_levels)  # at least one patch per level

    @data(
        *permute(
            dup(
                dict(
                    max_nbr_levels=3,
                    refinement="tagging",
                )
            ),
            expected_num_levels=3,
        ),
        *permute(dup(dict()), expected_num_levels=2),  # refinement boxes set later
    )
    @unpack
    def test_restarts(self, ndim, interp, simInput, expected_num_levels):
        print(f"test_restarts dim/interp:{ndim}/{interp}")

        simput = copy.deepcopy(simInput)

        for key in ["cells", "dl", "boundary_types"]:
            simput[key] = [simput[key]] * ndim

        if "refinement" not in simput:
            # three levels has issues with refinementboxes and possibly regridding
            b0 = [[10] * ndim, [19] * ndim]
            simput["refinement_boxes"] = {"L0": {"B0": b0}}
        else:  # https://github.com/LLNL/SAMRAI/issues/199
            # tagging can handle more than one timestep as it does not
            #  appear subject to regridding issues, so we make more timesteps
            #  to confirm simulations are still equivalent
            simput["time_step_nbr"] = 10

        # if restart time exists it "loads" from restart file
        #  otherwise just saves restart files based on timestamps
        assert "restart_time" not in simput["restart_options"]

        simput["interp_order"] = interp
        time_step = simput["time_step"]
        time_step_nbr = simput["time_step_nbr"]

        restart_idx = 4
        restart_time = time_step * restart_idx
        timestamps = [restart_time, time_step * time_step_nbr]

        # first simulation
        local_out = self.unique_diag_dir_for_test_case(f"{out}/test", ndim, interp)
        simput["restart_options"]["dir"] = local_out
        simput["restart_options"]["timestamps"] = [timestep * 4]
        simput["diag_options"]["options"]["dir"] = local_out
        ph.global_vars.sim = None
        ph.global_vars.sim = ph.Simulation(**simput)
        assert "restart_time" not in ph.global_vars.sim.restart_options
        model = setup_model()
        dump_all_diags(model.populations, timestamps=np.array(timestamps))
        Simulator(ph.global_vars.sim).run().reset()
        self.register_diag_dir_for_cleanup(local_out)
        diag_dir0 = local_out

        # second restarted simulation
        local_out = f"{local_out}_n2"
        simput["diag_options"]["options"]["dir"] = local_out
        simput["restart_options"]["restart_time"] = restart_time
        ph.global_vars.sim = None
        ph.global_vars.sim = ph.Simulation(**simput)
        assert "restart_time" in ph.global_vars.sim.restart_options
        model = setup_model()
        dump_all_diags(model.populations, timestamps=np.array(timestamps))
        Simulator(ph.global_vars.sim).run().reset()
        self.register_diag_dir_for_cleanup(local_out)
        diag_dir1 = local_out

        self.check_diags(
            diag_dir0, diag_dir1, model.populations, timestamps, expected_num_levels
        )

    @data(
        *permute(
            dup(
                dict(
                    max_nbr_levels=2,
                    refinement="tagging",
                )
            ),
            expected_num_levels=2,
        ),
    )
    @unpack
    def test_restarts_elapsed_time(self, ndim, interp, simInput, expected_num_levels):
        print(f"test_restarts_elapsed_time dim/interp:{ndim}/{interp}")

        simput = copy.deepcopy(simInput)
        simput["time_step_nbr"] = 2  # forcing restart after first advance

        for key in ["cells", "dl", "boundary_types"]:
            simput[key] = [simput[key]] * ndim

        if "refinement" not in simput:
            raise RuntimeError("No refinement box version due to regridding!")

        # if restart time exists it "loads" from restart file
        #  otherwise just saves restart files based on timestamps
        assert "restart_time" not in simput["restart_options"]

        simput["interp_order"] = interp
        time_step = simput["time_step"]
        time_step_nbr = simput["time_step_nbr"]
        timestamps = [time_step * time_step_nbr]

        # first simulation
        local_out = self.unique_diag_dir_for_test_case(
            f"{out}/elapsed_test", ndim, interp
        )
        diag_dir0 = local_out
        diag_dir1 = f"{local_out}_n2"

        seconds = 1  # dump on first advance always!
        simput["restart_options"]["elapsed_timestamps"] = [
            datetime.timedelta(seconds=seconds)
        ]
        simput["restart_options"]["dir"] = diag_dir0
        simput["diag_options"]["options"]["dir"] = diag_dir0
        ph.global_vars.sim = None
        ph.global_vars.sim = ph.Simulation(**simput)
        self.assertEqual(
            [seconds], ph.global_vars.sim.restart_options["elapsed_timestamps"]
        )

        assert "restart_time" not in ph.global_vars.sim.restart_options
        model = setup_model()
        dump_all_diags(model.populations, timestamps=np.array(timestamps))

        # autodump false to ignore possible init dump
        simulator = Simulator(ph.global_vars.sim, auto_dump=False).initialize()

        sleep(5)
        simulator.advance().dump()  # should trigger restart on "restart_idx" advance
        simulator.advance().dump()
        simulator.reset()
        self.register_diag_dir_for_cleanup(diag_dir0)

        # second restarted simulation
        simput["diag_options"]["options"]["dir"] = diag_dir1
        simput["restart_options"]["restart_time"] = time_step
        ph.global_vars.sim = None
        del simput["restart_options"]["elapsed_timestamps"]
        ph.global_vars.sim = ph.Simulation(**simput)
        assert "restart_time" in ph.global_vars.sim.restart_options
        model = setup_model()
        dump_all_diags(model.populations, timestamps=np.array(timestamps))
        Simulator(ph.global_vars.sim).run().reset()
        self.register_diag_dir_for_cleanup(diag_dir1)

        self.check_diags(
            diag_dir0, diag_dir1, model.populations, timestamps, expected_num_levels
        )

    def test_mode_conserve(self, ndim=1, interp=1, simput=dup(simArgs)):
        print(f"test_mode_conserve dim/interp:{ndim}/{interp}")

        for key in ["cells", "dl", "boundary_types"]:
            simput[key] = [simput[key]] * ndim

        # first simulation
        local_out = self.unique_diag_dir_for_test_case(f"{out}/conserve", ndim, interp)
        self.register_diag_dir_for_cleanup(local_out)

        simput["restart_options"]["dir"] = local_out
        simput["restart_options"]["timestamps"] = [timestep * 4]
        ph.global_vars.sim = ph.Simulation(**simput)
        self.assertEqual(len(ph.global_vars.sim.restart_options["timestamps"]), 1)
        self.assertEqual(ph.global_vars.sim.restart_options["timestamps"][0], 0.004)
        setup_model()
        Simulator(ph.global_vars.sim).run().reset()

        # second simulation (not restarted)
        ph.global_vars.sim = None
        simput["restart_options"]["mode"] = "conserve"
        ph.global_vars.sim = ph.Simulation(**simput)
        self.assertEqual(len(ph.global_vars.sim.restart_options["timestamps"]), 0)

    def test_input_validation_trailing_slash(self):
        if cpp.mpi_size() > 1:
            return  # no need to test in parallel

        simulation_args = dup()
        simulation_args["restart_options"]["dir"] += (
            simulation_args["restart_options"]["dir"] + "//"
        )
        sim = ph.Simulation(**simulation_args)
        setup_model()
        Simulator(sim).run().reset()
        ph.global_vars.sim = None

    @data(
        ([timedelta(hours=1), timedelta(hours=2)], True),
        ([timedelta(minutes=60), timedelta(minutes=120)], True),
        ([timedelta(minutes=60), timedelta(minutes=119)], False),
        ([timedelta(minutes=60), timedelta(minutes=59)], False),
    )
    @unpack
    def test_elapsed_timestamps_are_valid(self, elapsed_timestamps, valid):
        simput = dup(dict())
        simput["restart_options"]["elapsed_timestamps"] = elapsed_timestamps

        try:
            ph.global_vars.sim = None
            ph.Simulation(**simput.copy())
            self.assertTrue(valid)
        except Exception:
            self.assertTrue(not valid)

    def test_advanced_restarts_options(self):
        """
        Dim / interp / etc are not relevant here
        """
        ndim, interp = 1, 1
        print("test_advanced_restarts_options")

        simput = copy.deepcopy(
            dup(
                dict(
                    cells=10,
                    time_step_nbr=10,
                    max_nbr_levels=1,
                    refinement="tagging",
                )
            )
        )

        simput["interp_order"] = interp
        time_step = simput["time_step"]
        time_step_nbr = simput["time_step_nbr"]

        timestamps = time_step * np.arange(time_step_nbr + 1)
        local_out = self.unique_diag_dir_for_test_case(f"{out}/test", ndim, interp)
        simput["restart_options"]["dir"] = local_out
        simput["restart_options"]["keep_last"] = 3
        simput["restart_options"]["timestamps"] = timestamps

        ph.global_vars.sim = None
        ph.global_vars.sim = ph.Simulation(**simput)
        setup_model()
        Simulator(ph.global_vars.sim).run().reset()
        self.register_diag_dir_for_cleanup(local_out)

        simput["restart_options"]["restart_time"] = "auto"
        self.assertEqual(0.01, ph.restarts.restart_time(simput["restart_options"]))

        dirs = []
        for path_object in Path(local_out).iterdir():
            if path_object.is_dir():
                try:
                    dirs.append(float(path_object.name))
                except ValueError:
                    ...  # skip

        dirs = sorted(dirs)
        for i, idx in enumerate(range(8, 11)):
            self.assertAlmostEqual(dirs[i], time_step * idx)


if __name__ == "__main__":
    unittest.main()
