import numbers

import numpy as np

from .dict_populator import DictPopulator


def populateDict(sim):
    dp = DictPopulator()

    if sim.refinement == "tagging":
        dp.add_string("simulation/AMR/refinement/tagging/hybrid_method", "default")

    dp.add_string("simulation/algo/ion_updater/pusher/name", sim.particle_pusher)

    dp.add_double("simulation/algo/ohm/resistivity", sim.resistivity)
    dp.add_double("simulation/algo/ohm/hyper_resistivity", sim.hyper_resistivity)
    dp.add_enum_int("simulation/algo/ohm/hyper_mode", "HyperMode", sim.hyper_mode)

    _populate_populations(dp, sim)
    _populate_electromagnetic(dp, sim)
    _populate_electrons(dp, sim)


def _populate_populations(dp, sim):
    init_model = sim.model
    modelDict = init_model.model_dict

    if init_model.nbr_populations() < 0:
        raise RuntimeError("Number of populations cannot be negative")
    dp.add_size_t("simulation/ions/nbrPopulations", init_model.nbr_populations())

    partinit = "particle_initializer"
    for pop_index, pop in enumerate(init_model.populations):
        pop_path = "simulation/ions/pop"
        partinit_path = pop_path + f"{pop_index:d}/" + partinit + "/"
        d = modelDict[pop]
        dp.add_string(pop_path + f"{pop_index:d}/name", pop)
        dp.add_double(pop_path + f"{pop_index:d}/mass", d["mass"])
        dp.add_string(partinit_path + "name", "maxwellian")

        dp.add_init_function(sim.ndim, partinit_path + "density", d["density"])
        for k in ["x", "y", "z"]:
            dp.add_init_function(
                sim.ndim, partinit_path + f"bulk_velocity_{k}", d[f"v{k}"]
            )
            dp.add_init_function(
                sim.ndim, partinit_path + f"thermal_velocity_{k}", d[f"vth{k}"]
            )
        dp.add_double(partinit_path + "charge", d["charge"])
        dp.add_string(partinit_path + "basis", "cartesian")
        if "init" in d and "seed" in d["init"]:
            dp.add_optional_size_t(partinit_path + "init/seed", d["init"]["seed"])

        dp.add_int(partinit_path + "nbr_part_per_cell", d["nbrParticlesPerCell"])
        dp.add_double(partinit_path + "density_cut_off", d["density_cut_off"])


def _populate_electromagnetic(dp, sim):
    modelDict = sim.model.model_dict

    dp.add_string("simulation/electromag/name", "EM")
    dp.add_string("simulation/electromag/electric/name", "E")
    dp.add_string("simulation/electromag/magnetic/name", "B")

    maginit_path = "simulation/electromag/magnetic/initializer/"
    dp.add_init_function(sim.ndim, maginit_path + "x_component", modelDict["bx"])
    dp.add_init_function(sim.ndim, maginit_path + "y_component", modelDict["by"])
    dp.add_init_function(sim.ndim, maginit_path + "z_component", modelDict["bz"])


def _populate_electrons(dp, sim):
    if sim.electrons is None:
        raise RuntimeError("Error - no electrons registered to this Simulation")

    closure = sim.electrons.closure
    # vacuum is not supported by the polytropic closure; for the others the check is asked for
    min_density = sim.min_initial_density
    if min_density is None:
        min_density = 0.01 if closure.name() == "polytropic" else 0.0
    sim.model.validate_density(sim, min_density)

    entries = list(sim.electrons.dict_path())
    if closure.name() == "polytropic":
        entries.append(("electrons/pressure_closure/s", _initial_entropy(sim, closure)))

    for path, value in entries:
        if isinstance(value, str):
            dp.add_string("simulation/" + path, value)
        elif isinstance(value, numbers.Real) and not isinstance(value, bool):
            dp.add_double("simulation/" + path, value)
        elif callable(value):
            dp.add_init_function(sim.ndim, "simulation/" + path, value)
        else:
            raise TypeError(
                f"electrons entry '{path}' should be a str, int, float or callable,"
                f" got {type(value).__name__}"
            )


def _initial_entropy(sim, closure):
    """
    s = Pe/n^gamma from the user's Pe and ion densities: the initial entropy of
    the polytropic closure, free of the noise of the deposited density
    """
    model = sim.model.model_dict
    pops = sim.model.populations

    def entropy(*x):
        n = sum(model[pop]["charge"] * model[pop]["density"](*x) for pop in pops)
        pe = closure.Pe(*x)
        return np.broadcast_to(pe / n**closure.gamma, np.shape(x[0])).astype(np.float64)

    return entropy
