from pyphare.pharein import clearDict, fn_wrapper


def populateDict():
    import pybindlibs.dictator as pp

    from .global_vars import sim as simulation

    def add_int(path, val):
        pp.add_int(path, int(val))

    def add_double(path, val):
        pp.add_double(path, float(val))

    add_string = pp.add_string
    addInitFunction = getattr(pp, "addInitFunction{:d}".format(simulation.ndim) + "D")

    add_double("time_step", simulation.timestep)
    add_double("final_time", simulation.final_time)

    add_int("nbr_cells/x", simulation.cells[0])
    add_double("mesh_size/x", simulation.dl[0])
    add_double("origin/x", simulation.origin[0])

    if simulation.ndim > 1:
        add_int("nbr_cells/y", simulation.cells[1])
        add_double("mesh_size/y", simulation.dl[1])
        add_double("origin/y", simulation.origin[1])

        if simulation.ndim > 2:
            add_int("nbr_cells/z", simulation.cells[2])
            add_double("mesh_size/z", simulation.dl[2])
            add_double("origin/z", simulation.origin[2])

    add_double("fv_method/resistivity", simulation.eta)
    add_double("fv_method/hyper_resistivity", simulation.nu)
    add_double("fv_method/heat_capacity_ratio", simulation.gamma)
    add_double("to_primitive/heat_capacity_ratio", simulation.gamma)
    add_double("to_conservative/heat_capacity_ratio", simulation.gamma)

    add_string("state/name", "state")

    d = simulation.model.model_dict

    addInitFunction("state/density/initializer", fn_wrapper(d["density"]))
    addInitFunction("state/velocity/initializer/x_component", fn_wrapper(d["vx"]))
    addInitFunction("state/velocity/initializer/y_component", fn_wrapper(d["vy"]))
    addInitFunction("state/velocity/initializer/z_component", fn_wrapper(d["vz"]))
    addInitFunction("state/magnetic/initializer/x_component", fn_wrapper(d["bx"]))
    addInitFunction("state/magnetic/initializer/y_component", fn_wrapper(d["by"]))
    addInitFunction("state/magnetic/initializer/z_component", fn_wrapper(d["bz"]))
    addInitFunction("state/pressure/initializer", fn_wrapper(d["p"]))
