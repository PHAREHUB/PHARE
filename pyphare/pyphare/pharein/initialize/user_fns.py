from .general import add_double, add_string, add_int, add_size_t, fn_wrapper
import pybindlibs.dictator as pp


def populateDict(sim):
    populate_electromag(sim)
    populate_particles(sim)


def populate_electromag(sim):
    addInitFunction = getattr(pp, "addInitFunction{:d}".format(sim.ndim) + "D")
    modelDict = sim.model.model_dict
    maginit_path = "simulation/electromag/magnetic/initializer/"
    addInitFunction(maginit_path + "x_component", fn_wrapper(modelDict["bx"]))
    addInitFunction(maginit_path + "y_component", fn_wrapper(modelDict["by"]))
    addInitFunction(maginit_path + "z_component", fn_wrapper(modelDict["bz"]))


def populate_particles(sim):
    addInitFunction = getattr(pp, "addInitFunction{:d}".format(sim.ndim) + "D")

    init_model = sim.model
    modelDict = init_model.model_dict

    if init_model.nbr_populations() < 0:
        raise RuntimeError("Number of populations cannot be negative")
    add_size_t("simulation/ions/nbrPopulations", init_model.nbr_populations())

    partinit = "particle_initializer"
    for pop_index, pop in enumerate(init_model.populations):
        pop_path = "simulation/ions/pop"
        partinit_path = pop_path + "{:d}/".format(pop_index) + partinit + "/"
        d = modelDict[pop]

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
