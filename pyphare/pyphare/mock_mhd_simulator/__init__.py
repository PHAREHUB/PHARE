import numpy as np


def is_scalar(arg):
    return not isinstance(arg, (list, tuple)) and not is_nd_array(arg)


def is_nd_array(arg):
    return isinstance(arg, np.ndarray)


# converts scalars to array of expected size
# converts lists to arrays
class py_fn_wrapper:
    def __init__(self, fn):
        self.fn = fn

    def __call__(self, *xyz):
        args = [np.asarray(arg) for arg in xyz]
        ret = self.fn(*args)
        if isinstance(ret, list):
            ret = np.asarray(ret)
        if is_scalar(ret):
            ret = np.full(len(args[-1]), ret)
        return ret


# Wrap calls to user init functions to turn C++ vectors to ndarrays,
#  and returned ndarrays to C++ span
class fn_wrapper(py_fn_wrapper):
    def __init__(self, fn):
        super().__init__(fn)

    def __call__(self, *xyz):
        from pyphare.cpp import cpp_etc_lib

        # convert numpy array to C++ SubSpan
        # couples vector init functions to C++
        return cpp_etc_lib().makePyArrayWrapper(super().__call__(*xyz))


def clearDict():
    """
    dict may contain dangling references from a previous simulation unless cleared
    """
    import pybindlibs.dictator as pp

    pp.stop()


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

    add_double("godunov/resistivity", simulation.eta)
    add_double("godunov/hyper_resistivity", simulation.nu)
    add_double("godunov/heat_capacity_ratio", simulation.gamma)
    add_string("godunov/terms", simulation.terms)
    add_string("godunov/reconstruction", simulation.reconstruction)
    add_string("godunov/limiter", simulation.limiter)
    add_string("godunov/riemann", simulation.riemann)
    add_double("to_primitive/heat_capacity_ratio", simulation.gamma)
    add_double("to_conservative/heat_capacity_ratio", simulation.gamma)
    add_string("integrator", simulation.integrator)

    add_string("state/name", "state")
    add_string("state1/name", "state1")
    add_string("state2/name", "state2")

    d = simulation.model.model_dict

    addInitFunction("state/density/initializer", fn_wrapper(d["density"]))
    addInitFunction("state/velocity/initializer/x_component", fn_wrapper(d["vx"]))
    addInitFunction("state/velocity/initializer/y_component", fn_wrapper(d["vy"]))
    addInitFunction("state/velocity/initializer/z_component", fn_wrapper(d["vz"]))
    addInitFunction("state/magnetic/initializer/x_component", fn_wrapper(d["bx"]))
    addInitFunction("state/magnetic/initializer/y_component", fn_wrapper(d["by"]))
    addInitFunction("state/magnetic/initializer/z_component", fn_wrapper(d["bz"]))
    addInitFunction("state/pressure/initializer", fn_wrapper(d["p"]))

    addInitFunction("state1/density/initializer", fn_wrapper(d["density"]))
    addInitFunction("state1/velocity/initializer/x_component", fn_wrapper(d["vx"]))
    addInitFunction("state1/velocity/initializer/y_component", fn_wrapper(d["vy"]))
    addInitFunction("state1/velocity/initializer/z_component", fn_wrapper(d["vz"]))
    addInitFunction("state1/magnetic/initializer/x_component", fn_wrapper(d["bx"]))
    addInitFunction("state1/magnetic/initializer/y_component", fn_wrapper(d["by"]))
    addInitFunction("state1/magnetic/initializer/z_component", fn_wrapper(d["bz"]))
    addInitFunction("state1/pressure/initializer", fn_wrapper(d["p"]))

    addInitFunction("state2/density/initializer", fn_wrapper(d["density"]))
    addInitFunction("state2/velocity/initializer/x_component", fn_wrapper(d["vx"]))
    addInitFunction("state2/velocity/initializer/y_component", fn_wrapper(d["vy"]))
    addInitFunction("state2/velocity/initializer/z_component", fn_wrapper(d["vz"]))
    addInitFunction("state2/magnetic/initializer/x_component", fn_wrapper(d["bx"]))
    addInitFunction("state2/magnetic/initializer/y_component", fn_wrapper(d["by"]))
    addInitFunction("state2/magnetic/initializer/z_component", fn_wrapper(d["bz"]))
    addInitFunction("state2/pressure/initializer", fn_wrapper(d["p"]))
