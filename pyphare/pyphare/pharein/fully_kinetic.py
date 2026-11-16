"""
Setup checks of a FullyKineticModel simulation, in PHARE units (B0, n0, V_A, d_i, 1/Omega_ci):
the plasma frequency of a species is omega_p = (c/V_A) |q| sqrt(n/m), its cyclotron frequency
Omega_c = |q| B / m and its Debye length v_th / omega_p.

Errors are setups the explicit scheme cannot run (unstable, or without a solution for E),
warnings are setups it runs inaccurately, errors as well with Simulation(strict=True).
"""

import warnings

import numpy as np

accurate_plasma_period = 0.2
accurate_cyclotron_period = 0.3
max_dx_over_debye = 3.0
# the pusher is not relativistic
max_speed_over_light = 0.3
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
    directions = "xyz"[: sim.ndim]
    dt = sim.time_stepper.time_step
    c = sim.c_over_va
    dl = np.asarray(sim.dl)
    dx = dl.min()

    def sample(fn):
        return np.asarray(fn(*coords), dtype=float) * np.ones_like(coords[0])

    B = np.sqrt(sum(sample(model.model_dict[f"b{k}"]) ** 2 for k in "xyz"))
    pops = {name: model.model_dict[name] for name in model.populations}
    for name, pop in pops.items():
        if pop["charge"] == 0 or pop["mass"] <= 0:
            raise ValueError(
                f"Error: population {name} needs a non zero charge and a positive mass"
            )

    density = {name: sample(pop["density"]) for name, pop in pops.items()}
    vbulk = {
        name: {k: sample(pop[f"v{k}"]) for k in "xyz"} for name, pop in pops.items()
    }
    vth = {
        name: {k: sample(pop[f"vth{k}"]) for k in "xyz"} for name, pop in pops.items()
    }

    # E starts at zero, so the Gauss residual div E - c^2 rho stays at -c^2 rho(t=0): a charge
    # density present at t=0 is invisible to the fields, a frozen neutralizing background.
    # Until E^0 can be initialized (Poisson), refuse a local charge separation.
    charge_density = sum(pop["charge"] * density[name] for name, pop in pops.items())
    abs_charge_density = sum(
        abs(pop["charge"]) * density[name] for name, pop in pops.items()
    )
    separation = np.abs(charge_density).max() / abs_charge_density.max()
    if separation > net_charge_tolerance:
        raise ValueError(
            f"Error: charge density {separation:.3g} (relative): E^0 = 0, an initial charge"
            " density would act as a frozen background"
        )

    def plasma_frequency2(name):
        pop = pops[name]
        return c**2 * pop["charge"] ** 2 * density[name] / pop["mass"]

    # Yee + leapfrog in a cold plasma: sin^2(w dt / 2) = (w_p dt / 2)^2 + c^2 dt^2 sum k_grid^2,
    # with the total plasma frequency
    w_p = np.sqrt(sum(plasma_frequency2(name) for name in pops)).max()
    stability = (w_p * dt / 2) ** 2 + (c * dt) ** 2 * np.sum(1.0 / dl**2)
    if stability >= 1:
        raise ValueError(
            f"Error: (w_p dt / 2)^2 + c^2 dt^2 sum 1/dx^2 = {stability:.3f} >= 1,"
            " the leapfrog is unstable"
        )

    problems = []
    if w_p * dt > accurate_plasma_period:
        problems.append(f"total plasma frequency times dt = {w_p * dt:.3f}")

    nbr_cells = np.prod(sim.cells)
    for name, pop in pops.items():
        omega_c = abs(pop["charge"]) * B.max() / pop["mass"]
        if omega_c * dt > accurate_cyclotron_period:
            problems.append(f"{name} cyclotron frequency times dt = {omega_c * dt:.3f}")

        # the displacement guard stops the run if a particle crosses a cell in a step;
        # the largest of N gaussian samples is about sqrt(2 ln N) thermal velocities
        tail = np.sqrt(2 * np.log(nbr_cells * pop["nbrParticlesPerCell"])) + 1
        for d, k in enumerate(directions):
            bulk = np.abs(vbulk[name][k]).max()
            if bulk * dt >= dl[d]:
                raise ValueError(
                    f"Error: {name} bulk velocity times dt = {bulk * dt:.3g} >= d{k}"
                )
            fastest = (np.abs(vbulk[name][k]) + tail * vth[name][k]).max()
            if fastest * dt >= dl[d]:
                problems.append(
                    f"{name} thermal velocity: (|v_bulk| + {tail:.1f} v_th) dt ="
                    f" {fastest * dt:.3g} >= d{k}, tail particles will stop the run"
                )

        speed = max(np.abs(vbulk[name][k]).max() for k in "xyz") + 3 * max(
            vth[name][k].max() for k in "xyz"
        )
        if speed / c > max_speed_over_light:
            problems.append(
                f"{name} (|v_bulk| + 3 v_th) / c = {speed / c:.2f} > {max_speed_over_light},"
                " the pusher is not relativistic, compare with the speed of light"
            )

    # finite grid heating is driven by the lightest species
    lightest = min(pop["mass"] for pop in pops.values())
    inverse_debye2 = sum(
        plasma_frequency2(name) / vth[name]["x"] ** 2
        for name, pop in pops.items()
        if pop["mass"] == lightest
    )
    debye = 1 / np.sqrt(inverse_debye2.max())
    if dx > max_dx_over_debye * debye:
        problems.append(
            f"Debye length {debye:.3g} of the lightest species: dx / lambda_D ="
            f" {dx / debye:.1f} > {max_dx_over_debye:g}, expect numerical heating"
        )

    for problem in problems:
        if sim.strict:
            raise RuntimeError(f"FullyKineticModel: {problem}")
        warnings.warn(f"FullyKineticModel: {problem}", UserWarning, stacklevel=2)
