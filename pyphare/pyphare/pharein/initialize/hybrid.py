import itertools

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
    _populate_collisions(dp, sim)

    _populate_populations(dp, sim)
    _populate_electromagnetic(dp, sim)
    _populate_electrons(dp, sim)


def _populate_collisions(dp, sim):
    if sim.collisions is None:
        return
    path = "simulation/algo/collisions/"
    collisions = sim.collisions
    dp.add_double(path + "nu0", collisions["nu0"])
    dp.add_int(path + "every", collisions["every"])
    dp.add_enum_int(path + "scheme", "ScatteringLaw", collisions["scheme"])
    dp.add_bool(path + "exact_conservation", collisions["exact_conservation"])
    dp.add_size_t(path + "seed", collisions["seed"])
    _warn_on_collision_accuracy(sim)


def _warn_on_collision_accuracy(sim, threshold=0.02):
    """
    s at the thermal relative speed, on a sample of the initial profiles of each pair of
    populations: the relaxation rates are off by roughly 0.23 sqrt(s) (see collisions.hpp)
    """
    collisions = sim.collisions
    model = sim.model.model_dict
    coords = np.meshgrid(
        *[
            np.linspace(0, n * dl, min(n, 16), endpoint=False)
            for n, dl in zip(sim.cells, sim.dl)
        ],
        indexing="ij",
    )
    coords = [c.ravel() for c in coords]

    def thermal_speed_squared(pop):
        return sum(np.asarray(model[pop][f"vth{k}"](*coords)) ** 2 for k in "xyz") / 3.0

    worst, worst_pair = 0.0, None
    for a, b in itertools.combinations_with_replacement(sim.model.populations, 2):
        ma, mb = model[a]["mass"], model[b]["mass"]
        za, zb = model[a]["charge"], model[b]["charge"]
        density = np.maximum(model[a]["density"](*coords), model[b]["density"](*coords))
        u = np.sqrt(thermal_speed_squared(a) + thermal_speed_squared(b))
        s = (
            collisions["nu0"]
            * (za * zb) ** 2
            * density
            * sim.time_stepper.time_step
            * collisions["every"]
            / ((ma * mb / (ma + mb)) ** 2 * u**3)
        )
        if s.max() > worst:
            worst, worst_pair = s.max(), (a, b)

    if worst > threshold:
        print(
            f"WARNING: collisions: s at the thermal speed reaches {worst:.2g} for {worst_pair},"
            f" relaxation rates off by roughly {0.23 * np.sqrt(worst):.0%};"
            f" reduce the time step or 'every' to get below {threshold}"
        )


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

    for item in sim.electrons.dict_path():
        if isinstance(item[1], str):
            dp.add_string("simulation/" + item[0], item[1])
        else:
            dp.add_double("simulation/" + item[0], item[1])
