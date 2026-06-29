import pybindlibs.dictator as pp

from .general import add_double, add_int, add_string, fn_wrapper


def populateDict(sim):
    addInitFunction = getattr(pp, "addInitFunction{:d}".format(sim.ndim) + "D")

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
    # static background field B0 of the B = B0 + B1 split (defaults to zero)
    addInitFunction(
        "simulation/mhd_state/external_magnetic/initializer/x_component",
        fn_wrapper(modelDict["b0x"]),
    )
    addInitFunction(
        "simulation/mhd_state/external_magnetic/initializer/y_component",
        fn_wrapper(modelDict["b0y"]),
    )
    addInitFunction(
        "simulation/mhd_state/external_magnetic/initializer/z_component",
        fn_wrapper(modelDict["b0z"]),
    )
    addInitFunction(
        "simulation/mhd_state/pressure/initializer", fn_wrapper(modelDict["p"])
    )

    # vector-potential init (2D): B0 = curl(a0z z_hat), B1 = curl(a1z z_hat). The mode strings
    # tell the C++ side which fields to build from their potential instead of from the component
    # functions above. The potential functions are always provided (zero by default) so the dict
    # keys exist regardless of mode.
    add_string("simulation/mhd_state/b0_init_mode", modelDict["b0_init_mode"])
    add_string("simulation/mhd_state/b1_init_mode", modelDict["b1_init_mode"])
    addInitFunction(
        "simulation/mhd_state/external_magnetic/initializer/potential_z",
        fn_wrapper(modelDict["a0z"]),
    )
    addInitFunction(
        "simulation/mhd_state/perturbed_magnetic/initializer/potential_z",
        fn_wrapper(modelDict["a1z"]),
    )
