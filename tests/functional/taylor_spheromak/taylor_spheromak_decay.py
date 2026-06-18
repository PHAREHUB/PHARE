#
#
#


import numpy as np
import pyphare.pharein as ph


ph.NO_GUI()


cells = (101, 101, 101)
dl = (0.1, 0.1, 0.1)
dx, dy, dz = dl

diag_outputs = "phare_outputs/test/taylor_spheromak"
time_step = 0.001
final_time = 3

restart_options = {
    "dir": "checkpoints",
    "mode": "overwrite",
    "timestamps": [],
    "restart_time": "auto",
}
max_nbr_levels = 1

timestamps = np.arange(0, final_time + time_step, final_time / 100)


def get_taylor_spheromak(x, y, z, cx, cy, cz, R=1.0, B0=1.0):
    from scipy.special import spherical_jn

    # 1. Calculate relative coordinates and radial distance
    dx, dy, dz = x - cx, y - cy, z - cz
    r = np.sqrt(dx**2 + dy**2 + dz**2)

    # 2. Create the Mask (The "Neutral Zone" boundary)
    # This prevents the ValueError by handling the whole array at once
    mask = (r <= R) & (r > 0)

    # Initialize output arrays with zeros
    bx = np.zeros_like(x)
    by = np.zeros_like(y)
    bz = np.zeros_like(z)

    # 3. Calculate only for points inside the "Ball of Potential"
    k = 4.493 / R
    rho = k * r[mask]

    # Spherical Bessel j1(rho)
    j1_rho = spherical_jn(1, rho)

    # Force-Free projection factor
    # In a discrete substrate, this is the "Integer Gearing"
    common_factor = B0 * j1_rho / r[mask]

    # 4. Populate the vector components
    bx[mask] = common_factor * (dz[mask] * dx[mask] / r[mask] - dy[mask])
    by[mask] = common_factor * (dz[mask] * dy[mask] / r[mask] + dx[mask])
    bz[mask] = common_factor * (-(dx[mask] ** 2 + dy[mask] ** 2) / r[mask])

    return bx, by, bz


def b3(sim, x, y, z):
    L = sim.simulation_domain()
    mid = np.array(L) / 2
    cx = mid[0]
    cy = mid[1]
    cz = mid[2]
    bx, by, bz = get_taylor_spheromak(x, y, z, cx, cy, cz)
    return bx, by, bz


def config():
    sim = ph.Simulation(
        time_step=time_step,
        final_time=final_time,
        dl=dl,
        cells=cells,
        refinement="tagging",
        max_nbr_levels=max_nbr_levels,
        hyper_resistivity=0.001,
        resistivity=0.001,
        diag_options={
            "format": "pharevtkhdf",
            "options": {"dir": diag_outputs},
        },
        restart_options=restart_options,
        strict=False,
    )

    def density(x, y, z):
        return 0.5

    def b(x, y, z):
        return b3(sim, x, y, z)

    def bx(x, y, z):
        return b(x, y, z)[0]

    def by(x, y, z):
        return b(x, y, z)[1]

    def bz(x, y, z):
        return b(x, y, z)[2]

    def vxyz(x, y, z):
        return 0.0

    def vthxyz(x, y, z):
        return 0.00001

    C = "xyz"
    vvv = {
        **{f"vbulk{c}": vxyz for c in C},
        **{f"vth{c}": vthxyz for c in C},
        "nbr_part_per_cell": 100,
    }
    protons = {
        "charge": 1,
        "mass": 1,
        "density": density,
        **vvv,
        "init": {"seed": 12334},
    }
    ph.MaxwellianFluidModel(bx=bx, by=by, bz=bz, protons=protons)
    ph.ElectronModel(closure="isothermal", Te=0.0)

    for quantity in ["mass_density", "bulkVelocity"]:
        ph.FluidDiagnostics(quantity=quantity, write_timestamps=timestamps)

    ph.FluidDiagnostics(
        quantity="density", write_timestamps=timestamps, population_name="protons"
    )

    for quantity in ["E", "B"]:
        ph.ElectromagDiagnostics(quantity=quantity, write_timestamps=timestamps)

    return sim


if ph.PHARE_EXE:
    config()
elif __name__ == "__main__":
    from pyphare.simulator.simulator import Simulator, startMPI

    startMPI()
    Simulator(config()).run()
