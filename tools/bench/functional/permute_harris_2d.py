#
#
# how to use
"""
PHARE_SCOPE_TIMING=1 python3 -Ou tools/bench/functional/permute_harris_2d.py
python3 tools/python3/bench.py print_summary -i .phare_bench/20260528-151026.tar.gz > perms
python3 -c "from tools.bench.functional.permutor import sort_summaries as sort; sort('perms')"
python3 tools/bench/functional/plot_permutations.py
"""

import numpy as np
from copy import deepcopy

import pyphare.pharein as ph
from pyphare.simulator.simulator import startMPI

### test defaults
ndim = 2
dl = [0.4] * ndim
time_step = 0.001
final_time = 0.001
cells = [100, 100]
ppc = 100


def _cells():
    ret = []
    for by in [1]:
        ret.append(deepcopy(cells))
        ret[-1][0] *= by
    return ret


permutables = [
    ("interp_order", [1]),
    ("cells", _cells()),
    ("ppc", [100]),
    ("max_nbr_levels", [3]),
    ("tag_buffer", [3, 4, 5, 10]),  #
    ("tagging_threshold", [0.1, 0.2, 0.3, 0.4]),  #
    ("tile_size", [3, 4, 5, 6, 7, 8]),  #
]


def config(**kwargs):
    L = 0.5

    sim = ph.Simulation(
        refinement="tagging",
        interp_order=kwargs.get("interp_order", 1),
        time_step=kwargs.get("time_step", time_step),
        final_time=kwargs.get("final_time", final_time),
        dl=kwargs.get("dl", dl),
        cells=kwargs.get("cells", cells),
        max_nbr_levels=kwargs.get("max_nbr_levels", 1),
        tag_buffer=kwargs.get("tag_buffer", 3),
        nesting_buffer=kwargs.get("nesting_buffer", 1),
        tagging_threshold=kwargs.get("tagging_threshold", 0.1),
        hyper_resistivity=kwargs.get("hyper_resistivity", 0.008),
        hyper_mode=kwargs.get("hyper_mode", "spatial"),
        resistivity=kwargs.get("resistivity", 0.001),
        write_reports=False,
        strict=False,
        clustering={"method": "tile", "tile_size": kwargs.get("tile_size", 8)},
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
        "nbr_part_per_cell": 100,
    }

    ph.MaxwellianFluidModel(
        bx=bx,
        by=by,
        bz=bz,
        protons={"charge": 1, "density": density, **vvv, "init": {"seed": 12334}},
    )
    ph.ElectronModel(closure="isothermal", Te=0.0)

    return sim


if __name__ == "__main__":
    from permutor import execute

    startMPI()
    execute(config, permutables)
