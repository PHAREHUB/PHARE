"""
Global invariants of a periodic hybrid run, read from the live simulation.

For kinetic ions and isothermal electrons in a periodic box, without
resistivity or hyper-resistivity, the continuous equations conserve:

- the number of ions of each species, N_s = sum_p w_p
- the mean magnetic field: the curl of a periodic E sums to zero
- the total ion momentum, P = sum_s m_s sum_p w_p v_p
- the total energy E = K + W_B + U_e with
    K   = 1/2 sum_s m_s sum_p w_p v_p^2   (ions)
    W_B = 1/2 sum_cells B^2                (magnetic)
    U_e = Te sum_cells Ne ln Ne            (isothermal electrons)

U_e is the bookkeeping potential of the pressure force: with
E = -Ve x B - grad(Ne Te)/Ne, d(K + W_B)/dt = -int Ve.grad(Ne Te), and electron
continuity gives d/dt int Ne ln Ne = int Ve.grad Ne. This needs div J = 0.

The particle number and the mean B are exact in the discrete scheme too. The
momentum and the energy are not: the particle-to-grid deposit and the
grid-to-particle gather are not adjoint, and particle noise heats the plasma.

All quantities are per unit cell volume, so they compare with each other on
one level. Only level 0 is summed: a refined level overwrites the coarse
fields under it, which needs a composite sum this module does not do.
Serial runs only: there is no MPI reduction.
"""

import numpy as np

# 1 = primal along that direction, Yee layout
_PRIMAL = {"Bx": (1, 0, 0), "By": (0, 1, 0), "Bz": (0, 0, 1), "n": (1, 1, 1)}


def owned(patch_data, quantity, ndim):
    """values of one patch on its own cells, plus the lower primal node only,
    so that primal nodes shared by two patches are counted once"""
    g = patch_data.nGhosts
    n = [int(u - l + 1) for l, u in zip(patch_data.lower, patch_data.upper)]
    shape = [n[d] + 2 * g + _PRIMAL[quantity][d] for d in range(ndim)]
    values = np.asarray(patch_data.data).reshape(shape)
    return values[tuple(slice(g, g + n[d]) for d in range(ndim))]


def invariants(simulator, masses, Te, ndim):
    """sums over level 0. masses = {population name: mass}"""
    level = simulator.data_wrangler().getPatchLevel(0)

    B_sum = np.zeros(3)
    W_B = 0.0
    getters = (level.getBx, level.getBy, level.getBz)
    for i, (name, getter) in enumerate(zip(("Bx", "By", "Bz"), getters)):
        for patch in getter():
            values = owned(patch, name, ndim)
            B_sum[i] += values.sum()
            W_B += 0.5 * (values**2).sum()

    U_e = 0.0
    n_min = np.inf
    for patch in level.getDensity():  # charge density = electron density
        density = owned(patch, "n", ndim)
        n_min = min(n_min, density.min())
        U_e += Te * (density * np.log(density)).sum()

    N = {}
    K = 0.0
    P = np.zeros(3)
    P_noise2 = 0.0
    particles = level.getParticles()
    for name, mass in masses.items():
        N[name] = 0.0
        for patch in particles[name]["domain"]:
            weight = np.asarray(patch.data.weight)
            velocity = np.asarray(patch.data.v).reshape(-1, 3)
            N[name] += weight.sum()
            P += mass * (weight[:, None] * velocity).sum(axis=0)
            K += 0.5 * mass * (weight * (velocity**2).sum(axis=1)).sum()
            P_noise2 += (mass * weight) ** 2 @ (velocity**2).sum(axis=1)

    return {
        "N": N,
        "P": P,
        # random-walk level of P for uncorrelated particles: the scale of the
        # momentum error a noisy discrete scheme can make
        "P_noise": np.sqrt(P_noise2),
        "K": K,
        "W_B": W_B,
        "U_e": U_e,
        "E": K + W_B + U_e,
        "B_sum": B_sum,
        "n_min": n_min,
    }
