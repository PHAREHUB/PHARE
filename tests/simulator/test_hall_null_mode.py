import glob
import unittest

import numpy as np

import pyphare.pharein as ph
from pyphare.pharesee.hierarchy import hierarchy_from
from pyphare.simulator.simulator import Simulator, startMPI

ph.NO_GUI()

N = 32
DL = 0.1
L = N * DL
DT = 2e-4
NSTEPS = 50
AMP = 1e-4
MODES = ((5, 11), (9, 4))


def run(hall, diag_dir, bx, by, bz):
    ph.global_vars.sim = None
    sim = ph.Simulation(
        model_options=["MHDModel"],
        reconstruction="WENOZ",
        limiter="None",
        riemann="Rusanov",
        mhd_timestepper="TVDRK3",
        hall=hall,
        eta=0.0,
        nu=0.0,
        gamma=5.0 / 3.0,
        max_mhd_level=1,
        time_step=DT,
        time_step_nbr=NSTEPS,
        cells=(N, N),
        dl=(DL, DL),
        diag_options={
            "format": "phareh5",
            "options": {"dir": diag_dir, "mode": "overwrite"},
        },
    )

    def p(x, y):
        return 0.01 + 0.0 * x

    ph.MHDModel(density=one, vx=zero, vy=zero, vz=zero, bx=bx, by=by, bz=bz, p=p)
    final = [NSTEPS * DT]
    ph.ElectromagDiagnostics(quantity="B", write_timestamps=final)
    ph.MHDDiagnostics(quantity="Etot", write_timestamps=final)
    Simulator(sim).run().reset()
    ph.global_vars.sim = None


def zero(x, y):
    return 0.0 * x


def one(x, y):
    return 1.0 + 0.0 * x


def perturbation(x, y, modes):
    return AMP * sum(np.cos(2 * np.pi * (mx * x + my * y) / L) for mx, my in modes)


def level0_data(h5_filename):
    hier = hierarchy_from(h5_filename=h5_filename)
    t = hier.times()[-1]
    out = {}
    for patch in hier.level(0, t).patches:
        for name, pd in patch.patch_datas.items():
            out.setdefault(name, []).append((tuple(patch.box.lower), pd.dataset[:]))
    return {name: sorted(rows, key=lambda r: r[0]) for name, rows in out.items()}


def max_rel_diff(a, b):
    m = 0.0
    for name in a:
        for (la, da), (lb, db) in zip(a[name], b[name]):
            assert la == lb
            m = max(m, float(np.max(np.abs(da - db))))
    return m / AMP


def hall_on_off_ratios(name, bx, by, bz):
    dirs = {
        hall: f"phare_outputs/hall_null_mode/{name}_hall_{hall}"
        for hall in (True, False)
    }
    for hall, d in dirs.items():
        run(hall, d, bx, by, bz)

    def load(d, pattern):
        return level0_data(glob.glob(f"{d}/{pattern}")[0])

    b_ratio = max_rel_diff(load(dirs[True], "EM_B.h5"), load(dirs[False], "EM_B.h5"))
    e_ratio = max_rel_diff(
        load(dirs[True], "*Etot*.h5"), load(dirs[False], "*Etot*.h5")
    )
    print(f"{name}: B ratio {b_ratio:.3e}  Etot ratio {e_ratio:.3e}", flush=True)
    return b_ratio, e_ratio


class HallNullModeTest(unittest.TestCase):
    def test_hall_vanishes_at_first_order_for_b0_perp_k(self):
        def bz(x, y):
            return 1.0 + perturbation(x, y, MODES)

        b_ratio, e_ratio = hall_on_off_ratios("perp", zero, zero, bz)
        self.assertLess(b_ratio, 1e-2)
        self.assertLess(e_ratio, 1e-2)

    def test_hall_acts_at_first_order_for_b0_parallel_k(self):
        def bx(x, y):
            return 1.0 + 0.0 * x

        def by(x, y):
            return perturbation(x, y, ((5, 0),))

        b_ratio, _ = hall_on_off_ratios("parallel", bx, by, zero)
        self.assertGreater(b_ratio, 1e-1)


if __name__ == "__main__":
    startMPI()
    unittest.main()
