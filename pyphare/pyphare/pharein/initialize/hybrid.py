def populateDict(sim):
    from . import DictPopulator

    dp = DictPopulator()

    if sim.refinement == "tagging":
        dp.add_string("simulation/AMR/refinement/tagging/hybrid_method", "default")

    dp.add_string("simulation/algo/ion_updater/pusher/name", sim.particle_pusher)

    dp.add_double("simulation/algo/ohm/resistivity", sim.resistivity)
    dp.add_double("simulation/algo/ohm/hyper_resistivity", sim.hyper_resistivity)
    dp.add_string("simulation/algo/ohm/hyper_mode", sim.hyper_mode)

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
        partinit_path = pop_path + "{:d}/".format(pop_index) + partinit + "/"
        d = modelDict[pop]
        dp.add_string(pop_path + "{:d}/name".format(pop_index), pop)
        dp.add_double(pop_path + "{:d}/mass".format(pop_index), d["mass"])
        dp.add_string(partinit_path + "name", "maxwellian")

        dp.add_init_function(sim.ndim, partinit_path + "density", d["density"])
        dp.add_init_function(sim.ndim, partinit_path + "bulk_velocity_x", d["vx"])
        dp.add_init_function(sim.ndim, partinit_path + "bulk_velocity_y", d["vy"])
        dp.add_init_function(sim.ndim, partinit_path + "bulk_velocity_z", d["vz"])
        dp.add_init_function(sim.ndim, partinit_path + "thermal_velocity_x", d["vthx"])
        dp.add_init_function(sim.ndim, partinit_path + "thermal_velocity_y", d["vthy"])
        dp.add_init_function(sim.ndim, partinit_path + "thermal_velocity_z", d["vthz"])
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
