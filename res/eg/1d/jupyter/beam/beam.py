#
#
#


import numpy as np

import pyphare.pharein as ph


def config(**kwargs):
    sim = ph.Simulation(
        final_time=100,
        time_step=0.001,
        boundary_types="periodic",
        cells=165,
        dl=0.2,
        hyper_resistivity=0.01,
        diag_options={
            "format": "phareh5",
            "options": {"dir": kwargs["diagdir"], "mode": "overwrite"},
        },
    )

    def densityMain(x):
        return 1.0

    def densityBeam(x):
        return 0.01

    def bx(x):
        return 1.0

    def by(x):
        return 0.0

    def bz(x):
        return 0.0

    def vB(x):
        U0 = kwargs.get("U0", 0.0)
        return U0

    def v0(x):
        return 0.0

    def vth(x):
        Ti = kwargs.get("Ti", 0.0)
        return np.sqrt(Ti)

    vMain = {
        "vbulkx": v0,
        "vbulky": v0,
        "vbulkz": v0,
        "vthx": vth,
        "vthy": vth,
        "vthz": vth,
    }

    vBulk = {
        "vbulkx": vB,
        "vbulky": v0,
        "vbulkz": v0,
        "vthx": vth,
        "vthy": vth,
        "vthz": vth,
    }

    ph.MaxwellianFluidModel(
        bx=bx,
        by=by,
        bz=bz,
        main={"charge": 1, "density": densityMain, **vMain},
        beam={"charge": 1, "density": densityBeam, **vBulk},
    )

    ph.ElectronModel(closure="isothermal", Te=kwargs.get("Te", 0.0))

    timestamps = np.arange(0, sim.final_time, 0.1)

    for quantity in ["B", "E"]:
        ph.ElectromagDiagnostics(
            quantity=quantity,
            write_timestamps=timestamps,
        )

    for pop_name in ["main", "beam"]:
        ph.ParticleDiagnostics(
            quantity="domain",
            population_name=pop_name,
            write_timestamps=timestamps,
        )

    return sim


def main():
    import sys
    from pyphare.simulator.simulator import Simulator

    if len(sys.argv) != 5:
        raise RunTimerError('This code needs 4 paramaters, "run_name", U0, Te, Ti')

    diagdir = sys.argv[1]
    U0 = float(sys.argv[2])
    Te = float(sys.argv[3])
    Ti = float(sys.argv[4])

    Simulator(
        config(diagdir=diagdir, U0=U0, Te=Te, Ti=Ti), print_one_line=True
    ).run().reset()
    ph.global_vars.sim = None


if __name__ == "__main__":
    main()
