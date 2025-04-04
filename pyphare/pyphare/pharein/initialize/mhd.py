import pybindlibs.dictator as pp

from .general import add_double, add_int, add_string, fn_wrapper


def populateDict(sim):

    addInitFunction = getattr(pp, "addInitFunction{:d}".format(sim.ndim) + "D")

    add_int("simulation/AMR/max_mhd_level", sim.max_nbr_levels)

    if sim.refinement == "tagging":
        add_string("simulation/AMR/refinement/tagging/mhd_method", "default")

    add_double("simulation/algo/fv_method/resistivity", sim.eta)
    add_double("simulation/algo/fv_method/hyper_resistivity", sim.nu)
    add_double("simulation/algo/fv_method/heat_capacity_ratio", sim.gamma)
    add_double("simulation/algo/fv_euler/heat_capacity_ratio", sim.gamma)
    add_double("simulation/algo/to_primitive/heat_capacity_ratio", sim.gamma)
    add_double("simulation/algo/to_conservative/heat_capacity_ratio", sim.gamma)

    add_string("simulation/mhd_state/name", "mhd_state")

    add_double(
        "simulation/mhd_state/to_conservative_init/heat_capacity_ratio", sim.gamma
    )

    init_model = sim.mhd_model
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
    addInitFunction(
        "simulation/mhd_state/pressure/initializer", fn_wrapper(modelDict["p"])
    )
