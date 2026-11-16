"""
Setup checks of a FullyKineticModel simulation, in PHARE units (B0, n0, V_A, d_i, 1/Omega_ci):
the plasma frequency of a species is omega_p = (c/V_A) |q| sqrt(n/m), its cyclotron frequency
Omega_c = |q| B / m and its Debye length v_th / omega_p.

Errors are setups the explicit scheme cannot run (unstable, or without a solution for E),
warnings are setups it runs inaccurately.
"""

import warnings

import numpy as np

# leapfrog stability limit of omega_p dt, and the accuracy one
unstable_plasma_period = 2.0
accurate_plasma_period = 0.2
accurate_cyclotron_period = 0.3
max_dx_over_debye = 3.0
# a particle at this many thermal velocities must still move less than a cell per step
thermal_tail = 4.0
net_charge_tolerance = 1e-6


def check_fully_kinetic_simulation(**kwargs):
    """checks needing only the Simulation arguments, raised at its construction"""
    if kwargs["time_stepper"].mode != "constant":
        raise ValueError(
            "Error: FullyKineticModel needs a constant time step, E is staggered by dt/2"
        )

    dt = kwargs["time_stepper"].time_step
    dl = np.asarray(kwargs["dl"])
    c = kwargs["c_over_va"]
    light = c * dt * np.sqrt(np.sum(1.0 / dl**2))
    if light >= 1:
        raise ValueError(
            f"Error: FullyKineticModel light CFL c dt sqrt(sum 1/dx^2) = {light:.3f} >= 1"
        )

    if kwargs["refinement"] == "tagging" or kwargs["max_nbr_levels"] != 1:
        raise ValueError("Error: FullyKineticModel supports one level only for now")

    if kwargs["restart_options"] is not None:
        raise ValueError("Error: FullyKineticModel does not support restarts yet")


def _cell_centers(sim):
    axes = [(np.arange(n) + 0.5) * d for n, d in zip(sim.cells, sim.dl)]
    return np.meshgrid(*axes, indexing="ij")


def check_fully_kinetic(sim):
    """checks needing the populations, done before the C++ dictionary is filled"""
    if sim.load_balancer is not None:
        raise ValueError("Error: FullyKineticModel does not support load balancing yet")

    model = sim.model
    coords = _cell_centers(sim)
    dt = sim.time_stepper.time_step
    c = sim.c_over_va
    dx = min(sim.dl)

    def sample(fn):
        return np.asarray(fn(*coords), dtype=float) * np.ones_like(coords[0])

    B = np.sqrt(sum(sample(model.model_dict[f"b{k}"]) ** 2 for k in "xyz"))

    charge_density = np.zeros_like(coords[0])
    abs_charge_density = np.zeros_like(coords[0])
    problems = []

    for name in model.populations:
        pop = model.model_dict[name]
        q, m = pop["charge"], pop["mass"]
        if q == 0 or m <= 0:
            raise ValueError(
                f"Error: population {name} needs a non zero charge and a positive mass"
            )

        n = sample(pop["density"])
        charge_density += q * n
        abs_charge_density += abs(q) * n

        omega_p = c * abs(q) * np.sqrt(n.max() / m)
        if omega_p * dt >= unstable_plasma_period:
            raise ValueError(
                f"Error: {name} plasma frequency times dt = {omega_p * dt:.3f}"
                f" >= {unstable_plasma_period}, the leapfrog is unstable"
            )
        if omega_p * dt > accurate_plasma_period:
            problems.append(f"{name} plasma frequency times dt = {omega_p * dt:.3f}")

        omega_c = abs(q) * B.max() / m
        if omega_c * dt > accurate_cyclotron_period:
            problems.append(f"{name} cyclotron frequency times dt = {omega_c * dt:.3f}")

        vth = max(sample(pop[f"vth{k}"]).max() for k in "xyz")
        if thermal_tail * vth * dt >= dx:
            problems.append(
                f"{name} thermal velocity: {thermal_tail:g} v_th dt = "
                f"{thermal_tail * vth * dt:.3g} >= dx, tail particles will stop the run"
            )

        if q < 0:
            vth_min = min(sample(pop[f"vth{k}"]).min() for k in "xyz")
            debye = vth_min / omega_p
            if dx > max_dx_over_debye * debye:
                problems.append(
                    f"{name} Debye length {debye:.3g}: dx / lambda_D = {dx / debye:.1f}"
                    f" > {max_dx_over_debye:g}, expect numerical heating"
                )

    net = abs(charge_density.sum()) / abs_charge_density.sum()
    if net > net_charge_tolerance:
        raise ValueError(
            f"Error: net charge {net:.3g} (relative), Gauss's law has no periodic solution"
        )

    for problem in problems:
        warnings.warn(f"FullyKineticModel: {problem}", UserWarning, stacklevel=2)
