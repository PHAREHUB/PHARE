def populateDict(sim):
    from . import DictPopulator

    dp = DictPopulator()

    dp.add_int("simulation/AMR/max_mhd_level", sim.max_mhd_level)

    if sim.refinement == "tagging":
        dp.add_string("simulation/AMR/refinement/tagging/mhd_method", "default")

    _populate_algo(dp, sim)
    _populate_state(dp, sim)


def _populate_algo(dp, sim):
    dp.add_double("simulation/algo/fv_method/resistivity", sim.eta)
    dp.add_double("simulation/algo/fv_method/hyper_resistivity", sim.nu)
    dp.add_double("simulation/algo/fv_method/heat_capacity_ratio", sim.gamma)
    dp.add_string("simulation/algo/fv_method/hyper_mode", sim.hyper_mode)
    dp.add_double("simulation/algo/to_primitive/heat_capacity_ratio", sim.gamma)
    dp.add_double("simulation/algo/to_conservative/heat_capacity_ratio", sim.gamma)
    dp.add_double("simulation/algo/constrained_transport/resistivity", sim.eta)
    dp.add_double("simulation/algo/constrained_transport/hyper_resistivity", sim.nu)
    dp.add_string("simulation/algo/constrained_transport/hyper_mode", sim.hyper_mode)


def _populate_state(dp, sim):
    modelDict = sim.model.model_dict

    dp.add_string("simulation/mhd_state/name", "mhd_state")
    dp.add_double(
        "simulation/mhd_state/to_conservative_init/heat_capacity_ratio", sim.gamma
    )

    dp.add_init_function(
        sim.ndim, "simulation/mhd_state/density/initializer", modelDict["density"]
    )

    velinit_path = "simulation/mhd_state/velocity/initializer/"
    maginit_path = "simulation/mhd_state/magnetic/initializer/"
    for k in ["x", "y", "z"]:
        dp.add_init_function(
            sim.ndim, velinit_path + f"{k}_component", modelDict[f"v{k}"]
        )
        dp.add_init_function(
            sim.ndim, maginit_path + f"{k}_component", modelDict[f"b{k}"]
        )

    dp.add_init_function(
        sim.ndim, "simulation/mhd_state/pressure/initializer", modelDict["p"]
    )
