import os

import numpy as np
import pybindlibs.dictator as pp
from pyphare.core.phare_utilities import is_scalar
from pyphare.pharein.load_balancer import LoadBalancer
from pyphare.pharein.simulation import deserialize as deserialize_sim
from pyphare.pharein.simulation import serialize as serialize_sim


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


def populateDict(sim):

    add_string("simulation/name", "simulation_test")
    add_int("simulation/dimension", sim.ndim)

    if sim.smallest_patch_size is not None:
        add_vector_int("simulation/AMR/smallest_patch_size", sim.smallest_patch_size)
    if sim.largest_patch_size is not None:
        add_vector_int("simulation/AMR/largest_patch_size", sim.largest_patch_size)

    add_string("simulation/grid/layout_type", sim.layout)
    add_int("simulation/grid/nbr_cells/x", sim.cells[0])
    add_double("simulation/grid/meshsize/x", sim.dl[0])
    add_double("simulation/grid/origin/x", sim.origin[0])
    add_string("simulation/grid/boundary_type/x", sim.boundary_types[0])

    if sim.ndim > 1:
        add_int("simulation/grid/nbr_cells/y", sim.cells[1])
        add_double("simulation/grid/meshsize/y", sim.dl[1])
        add_double("simulation/grid/origin/y", sim.origin[1])
        add_string("simulation/grid/boundary_type/y", sim.boundary_types[1])

        if sim.ndim > 2:
            add_int("simulation/grid/nbr_cells/z", sim.cells[2])
            add_double("simulation/grid/meshsize/z", sim.dl[2])
            add_double("simulation/grid/origin/z", sim.origin[2])
            add_string("simulation/grid/boundary_type/z", sim.boundary_types[2])

    add_int("simulation/interp_order", sim.interp_order)
    add_int("simulation/refined_particle_nbr", sim.refined_particle_nbr)
    add_double("simulation/time_step", sim.time_step)
    add_int("simulation/time_step_nbr", sim.time_step_nbr)

    add_string("simulation/AMR/clustering", sim.clustering)
    add_vector_int("simulation/AMR/nesting_buffer", sim.nesting_buffer)
    add_int("simulation/AMR/tag_buffer", sim.tag_buffer)

    add_int("simulation/AMR/max_nbr_levels", sim.max_nbr_levels)

    add_int("simulation/AMR/max_mhd_level", sim.max_mhd_level)

    refinement_boxes = sim.refinement_boxes

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

    if refinement_boxes is not None and sim.refinement == "boxes":
        as_paths(refinement_boxes)
    elif sim.refinement == "tagging":
        add_string("simulation/AMR/refinement/tagging/method", "auto")
        # the two following params are hard-coded for now
        # they will become configurable when we have multi-models or several methods
        # per model
        add_double("simulation/AMR/refinement/tagging/threshold", sim.tagging_threshold)
    else:
        add_string(
            "simulation/AMR/refinement/tagging/method", "none"
        )  # integrator.h might want some looking at

    # load balancer block start
    lb = sim.load_balancer or LoadBalancer(active=False, _register=False)
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

    serialized_sim = serialize_sim(sim)

    #### adding diagnostics

    diag_path = "simulation/diagnostics/"
    for diag in list(sim.diagnostics.values()):
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
            if attr_key == "heat_capacity_ratio":
                add_double(
                    name_path + "/" + f"attribute_{attr_idx}_value",
                    diag.attributes[attr_key],
                )
            else:
                add_string(
                    name_path + "/" + f"attribute_{attr_idx}_value",
                    diag.attributes[attr_key],
                )

    if len(sim.diagnostics) > 0:
        if sim.diag_options is not None and "options" in sim.diag_options:
            add_string(diag_path + "filePath", sim.diag_options["options"]["dir"])
            if "mode" in sim.diag_options["options"]:
                add_string(diag_path + "mode", sim.diag_options["options"]["mode"])
            if "fine_dump_lvl_max" in sim.diag_options["options"]:
                add_int(
                    diag_path + "fine_dump_lvl_max",
                    sim.diag_options["options"]["fine_dump_lvl_max"],
                )
        else:
            add_string(diag_path + "filePath", "phare_output")
    #### diagnostics added

    #### adding restarts
    if sim.restart_options is not None:
        restart_options = sim.restart_options
        restarts_path = "simulation/restarts/"
        restart_file_path = "phare_outputs"

        if "dir" in restart_options:
            restart_file_path = restart_options["dir"]

        if "restart_time" in restart_options:
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
            if not sim.is_restartable_compared_to(deserialized_simulation):
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
