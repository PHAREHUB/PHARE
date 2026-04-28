#
#
#


import numpy as np

import pyphare.pharein as ph


time_step = 0.005
final_time = 100.0
dt = time_step * 100
timestamps = np.arange(0, final_time + dt, dt)


def density(x):
    return 1.0


def bx(x):
    return 1.0


def by(x):
    return 0.0


def bz(x):
    return 0.0


def T(x):
    return 0.125**2


def vWeak(x):
    from pyphare.pharein.global_vars import sim

    L = sim.simulation_domain()[0]
    x0 = 0.5 * L
    sigma = 2.0
    bubble = 0.08 * np.exp(-((x - x0) ** 2) / (2 * sigma**2))
    return bubble


def vNull(x):
    return 0.0


def vth(x):
    return np.sqrt(T(x))


vvv = {
    "vbulkx": vWeak,
    "vbulky": vNull,
    "vbulkz": vNull,
    "vthx": vth,
    "vthy": vth,
    "vthz": vth,
}


def config(**kwargs):
    sim = ph.Simulation(
        time_step=time_step,
        final_time=final_time,
        boundary_types="periodic",
        hyper_resistivity=0.001,
        cells=512,
        dl=0.25,
        diag_options={
            "format": "phareh5",
            "options": {"dir": kwargs["diagdir"], "mode": "overwrite"},
        },
    )

    ph.MaxwellianFluidModel(
        bx=bx,
        by=by,
        bz=bz,
        protons={"charge": 1, "density": density, "nbr_part_per_cell": 200, **vvv},
    )

    ph.ElectronModel(closure="isothermal", Te=kwargs.get("Te", 0.0))

    for quantity in ["E", "B"]:
        ph.ElectromagDiagnostics(quantity=quantity, write_timestamps=timestamps)

    for quantity in ["charge_density", "bulkVelocity"]:
        ph.FluidDiagnostics(quantity=quantity, write_timestamps=timestamps)

    for quantity in ["domain"]:
        ph.ParticleDiagnostics(
            quantity=quantity, write_timestamps=timestamps, population_name="protons"
        )
    return sim


if __name__ == "__main__":
    import sys
    from pyphare.simulator.simulator import Simulator

    if len(sys.argv) != 4:
        print('This code needs 3 paramaters, "run_name", Te, Ti')
    else:
        diagdir = sys.argv[1]
        Te = float(sys.argv[2])
        Ti = float(sys.argv[3])

    Simulator(config(diagdir=diagdir, Te=Te), print_one_line=True).run().reset()
    ph.global_vars.sim = None
