import pybindlibs.dictator as pp

from .general import add_bool, add_double, add_int, add_string, addSpaceTimeFunction, fn_wrapper


def populateDict(sim):
    addInitFunction = getattr(pp, "addInitFunction{:d}".format(sim.ndim) + "D")

    def addSTFunction(path, fn):
        addSpaceTimeFunction(path, fn, sim.ndim)

    add_int("simulation/AMR/max_mhd_level", sim.max_mhd_level)

    if sim.refinement == "tagging":
        add_string("simulation/AMR/refinement/tagging/mhd_method", "default")

    add_double("simulation/algo/fv_method/resistivity", sim.eta)
    add_double("simulation/algo/fv_method/hyper_resistivity", sim.nu)
    add_double("simulation/algo/fv_method/heat_capacity_ratio", sim.gamma)
    add_string("simulation/algo/fv_method/hyper_mode", sim.hyper_mode)
    add_double("simulation/algo/to_primitive/heat_capacity_ratio", sim.gamma)
    add_double("simulation/algo/to_conservative/heat_capacity_ratio", sim.gamma)
    add_double("simulation/algo/constrained_transport/resistivity", sim.eta)
    add_double("simulation/algo/constrained_transport/hyper_resistivity", sim.nu)
    add_string("simulation/algo/constrained_transport/hyper_mode", sim.hyper_mode)

    add_string("simulation/mhd_state/name", "mhd_state")

    add_double(
        "simulation/mhd_state/to_conservative_init/heat_capacity_ratio", sim.gamma
    )

    init_model = sim.model
    modelDict = init_model.model_dict

    addInitFunction(
        "simulation/mhd_state/density/initializer", fn_wrapper(modelDict["density"])
    )
    addInitFunction(
        "simulation/mhd_state/velocity/initializer/x_component",
        fn_wrapper(modelDict["vx"]),
    )
    addInitFunction(
        "simulation/mhd_state/velocity/initializer/y_component",
        fn_wrapper(modelDict["vy"]),
    )
    addInitFunction(
        "simulation/mhd_state/velocity/initializer/z_component",
        fn_wrapper(modelDict["vz"]),
    )
    addInitFunction(
        "simulation/mhd_state/magnetic/initializer/x_component",
        fn_wrapper(modelDict["bx"]),
    )
    addInitFunction(
        "simulation/mhd_state/magnetic/initializer/y_component",
        fn_wrapper(modelDict["by"]),
    )
    addInitFunction(
        "simulation/mhd_state/magnetic/initializer/z_component",
        fn_wrapper(modelDict["bz"]),
    )
    # background field B0 of the B = B0 + B1 split (defaults to zero). When time-dependent, the
    # component functions are space+time f(x,t) and are re-stamped each timestep by the C++ side.
    b0_time_dependent = modelDict["b0_time_dependent"]
    add_bool("simulation/mhd_state/external_magnetic/time_dependent", b0_time_dependent)
    b0_components_st = b0_time_dependent and modelDict["b0_init_mode"] == "components"
    addB0Function = addSTFunction if b0_components_st else (
        lambda path, fn: addInitFunction(path, fn_wrapper(fn))
    )
    addB0Function(
        "simulation/mhd_state/external_magnetic/initializer/x_component", modelDict["b0x"]
    )
    addB0Function(
        "simulation/mhd_state/external_magnetic/initializer/y_component", modelDict["b0y"]
    )
    addB0Function(
        "simulation/mhd_state/external_magnetic/initializer/z_component", modelDict["b0z"]
    )
    if b0_components_st:
        addSTFunction(
            "simulation/mhd_state/external_magnetic/derivative/x_component",
            modelDict["db0x_dt"],
        )
        addSTFunction(
            "simulation/mhd_state/external_magnetic/derivative/y_component",
            modelDict["db0y_dt"],
        )
        addSTFunction(
            "simulation/mhd_state/external_magnetic/derivative/z_component",
            modelDict["db0z_dt"],
        )
    addInitFunction(
        "simulation/mhd_state/pressure/initializer", fn_wrapper(modelDict["p"])
    )

    # vector-potential init: B0 = curl(a0), B1 = curl(a1) with full 3D vector potentials. The mode
    # strings tell the C++ side which fields to build from their potential instead of from the
    # component functions above. The potential components are always provided (zero by default) so
    # the dict keys exist regardless of mode; they are read as a VecFieldInitializer.
    add_string("simulation/mhd_state/b0_init_mode", modelDict["b0_init_mode"])
    add_string("simulation/mhd_state/b1_init_mode", modelDict["b1_init_mode"])
    # When B0 comes from a time-dependent potential, a0* are space+time f(x,t) and dB0/dt is built
    # on the C++ side as curl(dA/dt) from the derivative potential da0*_dt (same discrete curl).
    a0_st = b0_time_dependent and modelDict["b0_init_mode"] == "potential"
    addA0Function = addSTFunction if a0_st else (
        lambda path, fn: addInitFunction(path, fn_wrapper(fn))
    )
    addA0Function(
        "simulation/mhd_state/external_magnetic/initializer/potential/x_component",
        modelDict["a0x"],
    )
    addA0Function(
        "simulation/mhd_state/external_magnetic/initializer/potential/y_component",
        modelDict["a0y"],
    )
    addA0Function(
        "simulation/mhd_state/external_magnetic/initializer/potential/z_component",
        modelDict["a0z"],
    )
    if a0_st:
        addSTFunction(
            "simulation/mhd_state/external_magnetic/derivative/potential/x_component",
            modelDict["da0x_dt"],
        )
        addSTFunction(
            "simulation/mhd_state/external_magnetic/derivative/potential/y_component",
            modelDict["da0y_dt"],
        )
        addSTFunction(
            "simulation/mhd_state/external_magnetic/derivative/potential/z_component",
            modelDict["da0z_dt"],
        )
    addInitFunction(
        "simulation/mhd_state/perturbed_magnetic/initializer/potential/x_component",
        fn_wrapper(modelDict["a1x"]),
    )
    addInitFunction(
        "simulation/mhd_state/perturbed_magnetic/initializer/potential/y_component",
        fn_wrapper(modelDict["a1y"]),
    )
    addInitFunction(
        "simulation/mhd_state/perturbed_magnetic/initializer/potential/z_component",
        fn_wrapper(modelDict["a1z"]),
    )
