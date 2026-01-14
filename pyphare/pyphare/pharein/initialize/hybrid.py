import pybindlibs.dictator as pp

from .general import add_double, add_int, add_size_t, add_string, fn_wrapper


def populateDict(sim):

    addInitFunction = getattr(pp, "addInitFunction{:d}".format(sim.ndim) + "D")

    if sim.refinement == "tagging":
        add_string("simulation/AMR/refinement/tagging/hybrid_method", "default")

    add_string("simulation/algo/ion_updater/pusher/name", sim.particle_pusher)

    add_double("simulation/algo/ohm/resistivity", sim.resistivity)
    add_double("simulation/algo/ohm/hyper_resistivity", sim.hyper_resistivity)
    add_string("simulation/algo/ohm/hyper_mode", sim.hyper_mode)

    init_model = sim.maxwellian_fluid_model
    modelDict = init_model.model_dict

    if init_model.nbr_populations() < 0:
        raise RuntimeError("Number of populations cannot be negative")
    add_size_t("simulation/ions/nbrPopulations", init_model.nbr_populations())

    partinit = "particle_initializer"
    for pop_index, pop in enumerate(init_model.populations):
        pop_path = "simulation/ions/pop"
        partinit_path = pop_path + "{:d}/".format(pop_index) + partinit + "/"
        d = modelDict[pop]
        add_string(pop_path + "{:d}/name".format(pop_index), pop)
        add_double(pop_path + "{:d}/mass".format(pop_index), d["mass"])
        add_string(partinit_path + "name", "maxwellian")

        addInitFunction(partinit_path + "density", fn_wrapper(d["density"]))
        addInitFunction(partinit_path + "bulk_velocity_x", fn_wrapper(d["vx"]))
        addInitFunction(partinit_path + "bulk_velocity_y", fn_wrapper(d["vy"]))
        addInitFunction(partinit_path + "bulk_velocity_z", fn_wrapper(d["vz"]))
        addInitFunction(partinit_path + "thermal_velocity_x", fn_wrapper(d["vthx"]))
        addInitFunction(partinit_path + "thermal_velocity_y", fn_wrapper(d["vthy"]))
        addInitFunction(partinit_path + "thermal_velocity_z", fn_wrapper(d["vthz"]))
        add_double(partinit_path + "charge", d["charge"])
        add_string(partinit_path + "basis", "cartesian")
        if "init" in d and "seed" in d["init"]:
            pp.add_optional_size_t(partinit_path + "init/seed", d["init"]["seed"])

        add_int(partinit_path + "nbr_part_per_cell", d["nbrParticlesPerCell"])
        add_double(partinit_path + "density_cut_off", d["density_cut_off"])

    add_string("simulation/electromag/name", "EM")
    add_string("simulation/electromag/electric/name", "E")
    add_string("simulation/electromag/magnetic/name", "B")

    maginit_path = "simulation/electromag/magnetic/initializer/"
    addInitFunction(maginit_path + "x_component", fn_wrapper(modelDict["bx"]))
    addInitFunction(maginit_path + "y_component", fn_wrapper(modelDict["by"]))
    addInitFunction(maginit_path + "z_component", fn_wrapper(modelDict["bz"]))

    #### adding electrons
    if sim.electrons is None:
        raise RuntimeError("Error - no electrons registered to this Simulation")
    else:
        for item in sim.electrons.dict_path():
            if isinstance(item[1], str):
                add_string("simulation/" + item[0], item[1])
            else:
                add_double("simulation/" + item[0], item[1])
