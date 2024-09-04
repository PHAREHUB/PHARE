from .general import add_string, add_int

# import pybindlibs.dictator as pp


def populateDict(sim):
    init_model = sim.model
    partinit = "particle_initializer"
    for pop_index, pop in enumerate(init_model.populations):
        pop_path = "simulation/ions/pop"
        partinit_path = pop_path + "{:d}/".format(pop_index) + partinit + "/"

        add_string(partinit_path + "name", "samraih5")
        add_string(partinit_path + "filepath", sim.init_options["dir"])
        add_int(partinit_path + "mpi_size", sim.init_options.get("mpi_size", 1))
        add_int(partinit_path + "index", sim.init_options.get("index", 0))
