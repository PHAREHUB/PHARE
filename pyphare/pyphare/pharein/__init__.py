import os
import sys
import subprocess
import numpy as np

from pyphare.core.phare_utilities import is_scalar
from .uniform_model import UniformModel
from .maxwellian_fluid_model import MaxwellianFluidModel
from .electron_model import ElectronModel
from .diagnostics import (
    FluidDiagnostics,
    ElectromagDiagnostics,
    ParticleDiagnostics,
    MetaDiagnostics,
    InfoDiagnostics,
)
from .simulation import (
    Simulation,
    serialize as serialize_sim,
    deserialize as deserialize_sim,
)
from .load_balancer import LoadBalancer

__all__ = [
    "UniformModel",
    "MaxwellianFluidModel",
    "ElectronModel",
    "FluidDiagnostics",
    "ElectromagDiagnostics",
    "ParticleDiagnostics",
    "MetaDiagnostics",
    "InfoDiagnostics",
    "Simulation",
]

# This exists to allow a condition variable for when we are running PHARE from C++ via phare-exe
#  It is configured to "True" in pyphare/pyphare/pharein/init.py::get_user_inputs(jobname)
#    which is called from Pybind when phare-exe (or similar) is in use
PHARE_EXE = False

venv_path = os.environ.get("VIRTUAL_ENV")

if venv_path is not None:
    pythonexe = os.path.join(venv_path, "bin/python3")
    arg = "import sys; print(sys.path)"
    p = subprocess.run([pythonexe, "-c", arg], stdout=subprocess.PIPE)
    s = p.stdout.decode()
    s = (
        s.replace("[", "")
        .replace("]", "")
        .replace('"', "")
        .replace("\n", "")
        .replace("'", "")
    )
    s = s.split(",")[1:]
    pythonpath = [ss.strip() for ss in s]
    sys.path = sys.path + pythonpath


def NO_GUI():
    """prevents issues when command line only and no desktop etc"""
    import matplotlib as mpl

    mpl.use("Agg")


def getSimulation():
    from .global_vars import sim

    return sim


def _patch_data_ids(restart_file_dir):
    """
    for restarts we save samrai patch data ids to the restart files, which we access from here
    to tell samrai which patch datas to load from the restart file on restart
    """
    from pyphare.cpp import cpp_etc_lib

    return cpp_etc_lib().patch_data_ids(restart_file_dir)


def _serialized_simulation_string(restart_file_dir):
    from pyphare.cpp import cpp_etc_lib

    return cpp_etc_lib().serialized_simulation_string(restart_file_dir)


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
    from .global_vars import sim as simulation
    import pybindlibs.dictator as pp

    # pybind complains if receiving wrong type
    def add_int(path, val):
        pp.add_int(path, int(val))

    def add_bool(path, val):
        pp.add_bool(path, bool(val))

    def add_double(path, val):
        pp.add_double(path, float(val))

    def add_size_t(path, val):
        casted = int(val)
        if casted < 0:
            raise RuntimeError("pyphare.__init__::add_size_t received negative value")
        pp.add_size_t(path, casted)

    def add_vector_int(path, val):
        pp.add_vector_int(path, list(val))

    add_string = pp.add_string
    addInitFunction = getattr(pp, "addInitFunction{:d}".format(simulation.ndim) + "D")

    add_string("simulation/name", "simulation_test")
    add_int("simulation/dimension", simulation.ndim)

    if simulation.smallest_patch_size is not None:
        add_vector_int(
            "simulation/AMR/smallest_patch_size", simulation.smallest_patch_size
        )
    if simulation.largest_patch_size is not None:
        add_vector_int(
            "simulation/AMR/largest_patch_size", simulation.largest_patch_size
        )

    add_string("simulation/grid/layout_type", simulation.layout)
    add_int("simulation/grid/nbr_cells/x", simulation.cells[0])
    add_double("simulation/grid/meshsize/x", simulation.dl[0])
    add_string("simulation/grid/boundary_type/x", simulation.boundary_types[0])

    if simulation.ndim > 1:
        add_int("simulation/grid/nbr_cells/y", simulation.cells[1])
        add_double("simulation/grid/meshsize/y", simulation.dl[1])
        add_string("simulation/grid/boundary_type/y", simulation.boundary_types[1])

        if simulation.ndim > 2:
            add_int("simulation/grid/nbr_cells/z", simulation.cells[2])
            add_double("simulation/grid/meshsize/z", simulation.dl[2])
            add_string("simulation/grid/boundary_type/z", simulation.boundary_types[2])

    add_int("simulation/interp_order", simulation.interp_order)
    add_int("simulation/refined_particle_nbr", simulation.refined_particle_nbr)
    add_double("simulation/time_step", simulation.time_step)
    add_int("simulation/time_step_nbr", simulation.time_step_nbr)

    add_string("simulation/AMR/clustering", simulation.clustering)
    add_int("simulation/AMR/max_nbr_levels", simulation.max_nbr_levels)
    add_vector_int("simulation/AMR/nesting_buffer", simulation.nesting_buffer)

    add_int("simulation/AMR/tag_buffer", simulation.tag_buffer)

    refinement_boxes = simulation.refinement_boxes

    def as_paths(rb):
        add_int("simulation/AMR/refinement/boxes/nbr_levels/", len(rb.keys()))
        for level, boxes in rb.items():
            level_path = "simulation/AMR/refinement/boxes/" + level + "/"
            add_int(level_path + "nbr_boxes/", int(len(boxes)))
            for box_i, box in enumerate(boxes):
                box_id = "B" + str(box_i)
                lower = box.lower
                upper = box.upper
                box_lower_path_x = box_id + "/lower/x/"
                box_upper_path_x = box_id + "/upper/x/"
                add_int(level_path + box_lower_path_x, lower[0])
                add_int(level_path + box_upper_path_x, upper[0])
                if len(lower) >= 2:
                    box_lower_path_y = box_id + "/lower/y/"
                    box_upper_path_y = box_id + "/upper/y/"
                    add_int(level_path + box_lower_path_y, lower[1])
                    add_int(level_path + box_upper_path_y, upper[1])
                    if len(lower) == 3:
                        box_lower_path_z = box_id + "/lower/z/"
                        box_upper_path_z = box_id + "/upper/z/"
                        add_int(level_path + box_lower_path_z, lower[2])
                        add_int(level_path + box_upper_path_z, upper[2])

    if refinement_boxes is not None and simulation.refinement == "boxes":
        as_paths(refinement_boxes)
    elif simulation.refinement == "tagging":
        add_string("simulation/AMR/refinement/tagging/method", "auto")
        # the two following params are hard-coded for now
        # they will become configurable when we have multi-models or several methods
        # per model
        add_string("simulation/AMR/refinement/tagging/model", "HybridModel")
        add_string("simulation/AMR/refinement/tagging/method", "default")
        add_double(
            "simulation/AMR/refinement/tagging/threshold", simulation.tagging_threshold
        )
    else:
        add_string(
            "simulation/AMR/refinement/tagging/method", "none"
        )  # integrator.h might want some looking at

    add_string("simulation/algo/ion_updater/pusher/name", simulation.particle_pusher)

    add_double("simulation/algo/ohm/resistivity", simulation.resistivity)
    add_double("simulation/algo/ohm/hyper_resistivity", simulation.hyper_resistivity)
    add_string("simulation/algo/ohm/hyper_mode", simulation.hyper_mode)

    # load balancer block start
    lb = simulation.load_balancer or LoadBalancer(active=False, _register=False)
    base = "simulation/AMR/loadbalancing"
    add_bool(f"{base}/active", lb.active)
    add_string(f"{base}/mode", lb.mode)
    add_double(f"{base}/tolerance", lb.tol)

    # if mode==nppc, imbalance allowed
    add_bool(f"{base}/auto", lb.auto)
    add_size_t(f"{base}/next_rebalance", lb.next_rebalance)
    add_size_t(f"{base}/max_next_rebalance", lb.max_next_rebalance)
    add_size_t(
        f"{base}/next_rebalance_backoff_multiplier",
        lb.next_rebalance_backoff_multiplier,
    )

    # cadence based values
    add_size_t(f"{base}/every", lb.every)
    add_bool(f"{base}/on_init", lb.on_init)
    # load balancer block end

    init_model = simulation.model
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

    serialized_sim = serialize_sim(simulation)

    #### adding diagnostics

    diag_path = "simulation/diagnostics/"
    for diag in list(simulation.diagnostics.values()):
        diag.attributes["serialized_simulation"] = serialized_sim

        type_path = diag_path + diag.type + "/"
        name_path = type_path + diag.name
        add_string(name_path + "/" + "type", diag.type)
        add_string(name_path + "/" + "quantity", diag.quantity)
        add_size_t(name_path + "/" + "flush_every", diag.flush_every)
        pp.add_array_as_vector(
            name_path + "/" + "write_timestamps", diag.write_timestamps
        )
        pp.add_array_as_vector(
            name_path + "/" + "compute_timestamps", diag.compute_timestamps
        )

        add_size_t(name_path + "/" + "n_attributes", len(diag.attributes))
        for attr_idx, attr_key in enumerate(diag.attributes):
            add_string(name_path + "/" + f"attribute_{attr_idx}_key", attr_key)
            add_string(
                name_path + "/" + f"attribute_{attr_idx}_value",
                diag.attributes[attr_key],
            )

    if len(simulation.diagnostics) > 0:
        if simulation.diag_options is not None and "format" in simulation.diag_options:
            add_string(diag_path + "format", simulation.diag_options["format"])

        if simulation.diag_options is not None and "options" in simulation.diag_options:
            add_string(
                diag_path + "filePath", simulation.diag_options["options"]["dir"]
            )
            if "mode" in simulation.diag_options["options"]:
                add_string(
                    diag_path + "mode", simulation.diag_options["options"]["mode"]
                )
            if "fine_dump_lvl_max" in simulation.diag_options["options"]:
                add_int(
                    diag_path + "fine_dump_lvl_max",
                    simulation.diag_options["options"]["fine_dump_lvl_max"],
                )
        else:
            add_string(diag_path + "filePath", "phare_output")
    #### diagnostics added

    #### adding restarts
    if simulation.restart_options is not None:
        restart_options = simulation.restart_options
        restarts_path = "simulation/restarts/"
        restart_file_path = "phare_outputs"

        if "dir" in restart_options:
            restart_file_path = restart_options["dir"]

        if "restart_time" in restart_options and restart_options["restart_time"] > 0:
            from pyphare.cpp import cpp_etc_lib

            restart_time = restart_options["restart_time"]
            restart_file_load_path = cpp_etc_lib().restart_path_for_time(
                restart_file_path, restart_time
            )

            if not os.path.exists(restart_file_load_path):
                raise ValueError(
                    f"PHARE restart file not found for time {restart_time}"
                )

            deserialized_simulation = deserialize_sim(
                _serialized_simulation_string(restart_file_load_path)
            )
            if not simulation.is_restartable_compared_to(deserialized_simulation):
                raise ValueError(
                    "deserialized Restart simulation is incompatible with configured simulation parameters"
                )

            add_vector_int(
                restarts_path + "restart_ids", _patch_data_ids(restart_file_load_path)
            )
            add_string(restarts_path + "loadPath", restart_file_load_path)
            add_double(restarts_path + "restart_time", restart_time)

        if "mode" in restart_options:
            add_string(restarts_path + "mode", restart_options["mode"])

        add_string(restarts_path + "filePath", restart_file_path)

        if "elapsed_timestamps" in restart_options:
            pp.add_array_as_vector(
                restarts_path + "elapsed_timestamps",
                restart_options["elapsed_timestamps"],
            )

        if "timestamps" in restart_options:
            pp.add_array_as_vector(
                restarts_path + "write_timestamps", restart_options["timestamps"]
            )

        add_string(restarts_path + "serialized_simulation", serialized_sim)
    #### restarts added

    #### adding electrons
    if simulation.electrons is None:
        raise RuntimeError("Error - no electrons registered to this Simulation")
    else:
        for item in simulation.electrons.dict_path():
            if isinstance(item[1], str):
                add_string("simulation/" + item[0], item[1])
            else:
                add_double("simulation/" + item[0], item[1])
