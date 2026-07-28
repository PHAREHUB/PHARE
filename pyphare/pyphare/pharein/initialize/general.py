import os

import numpy as np
from pyphare.core.phare_utilities import is_scalar
from pyphare.pharein.load_balancer import LoadBalancer
from pyphare.pharein.simulation import deserialize as deserialize_sim
from pyphare.pharein.simulation import serialize as serialize_sim


def populateDict(sim):
    from . import DictPopulator

    dp = DictPopulator()

    populate_grid(dp, sim)
    populate_amr(dp, sim)
    populate_load_balancer(dp, sim)

    serialized_sim = serialize_sim(sim)
    populate_diagnostics(dp, sim, serialized_sim)
    populate_restarts(dp, sim, serialized_sim)


def populate_grid(dp, sim):
    dp.add_string("simulation/name", "simulation_test")
    dp.add_int("simulation/dimension", sim.ndim)

    dp.add_string("simulation/grid/layout_type", sim.layout)
    dp.add_int("simulation/grid/nbr_cells/x", sim.cells[0])
    dp.add_double("simulation/grid/meshsize/x", sim.dl[0])
    dp.add_string("simulation/grid/boundary_type/x", sim.boundary_types[0])

    if sim.ndim > 1:
        dp.add_int("simulation/grid/nbr_cells/y", sim.cells[1])
        dp.add_double("simulation/grid/meshsize/y", sim.dl[1])
        dp.add_string("simulation/grid/boundary_type/y", sim.boundary_types[1])

        if sim.ndim > 2:
            dp.add_int("simulation/grid/nbr_cells/z", sim.cells[2])
            dp.add_double("simulation/grid/meshsize/z", sim.dl[2])
            dp.add_string("simulation/grid/boundary_type/z", sim.boundary_types[2])

    dp.add_int("simulation/interp_order", sim.interp_order)
    dp.add_int("simulation/refined_particle_nbr", sim.refined_particle_nbr)
    dp.add_double("simulation/time_step", sim.time_step)
    dp.add_int("simulation/time_step_nbr", sim.time_step_nbr)
    dp.add_double("simulation/final_time", sim.final_time)


def populate_amr(dp, sim):
    if sim.smallest_patch_size is not None:
        dp.add_vector_int("simulation/AMR/smallest_patch_size", sim.smallest_patch_size)
    if sim.largest_patch_size is not None:
        dp.add_vector_int("simulation/AMR/largest_patch_size", sim.largest_patch_size)

    dp.add_string("simulation/AMR/clustering", sim.clustering)
    dp.add_vector_int("simulation/AMR/nesting_buffer", sim.nesting_buffer)
    dp.add_int("simulation/AMR/tag_buffer", sim.tag_buffer)

    dp.add_int("simulation/AMR/max_nbr_levels", sim.max_nbr_levels)

    dp.add_int("simulation/AMR/max_mhd_level", sim.max_mhd_level)

    refinement_boxes = sim.refinement_boxes

    if refinement_boxes is not None and sim.refinement == "boxes":
        _populate_refinement_boxes(dp, refinement_boxes)
    elif sim.refinement == "tagging":
        dp.add_string("simulation/AMR/refinement/tagging/method", "auto")
        # the two following params are hard-coded for now
        # they will become configurable when we have multi-models or several methods
        # per model
        dp.add_double(
            "simulation/AMR/refinement/tagging/threshold", sim.tagging_threshold
        )
    else:
        dp.add_string(
            "simulation/AMR/refinement/tagging/method", "none"
        )  # integrator.h might want some looking at


def populate_load_balancer(dp, sim):
    lb = sim.load_balancer or LoadBalancer(active=False, _register=False)
    base = "simulation/AMR/loadbalancing"
    dp.add_bool(f"{base}/active", lb.active)
    dp.add_string(f"{base}/mode", lb.mode)
    dp.add_double(f"{base}/tolerance", lb.tol)

    # if mode==nppc, imbalance allowed
    dp.add_bool(f"{base}/auto", lb.auto)
    dp.add_size_t(f"{base}/next_rebalance", lb.next_rebalance)
    dp.add_size_t(f"{base}/max_next_rebalance", lb.max_next_rebalance)
    dp.add_size_t(
        f"{base}/next_rebalance_backoff_multiplier",
        lb.next_rebalance_backoff_multiplier,
    )

    # cadence based values
    dp.add_size_t(f"{base}/every", lb.every)
    dp.add_bool(f"{base}/on_init", lb.on_init)


def populate_diagnostics(dp, sim, serialized_sim):
    diag_path = "simulation/diagnostics/"
    for diag in list(sim.diagnostics.values()):
        diag.attributes["serialized_simulation"] = serialized_sim

        name_path = diag_path + diag.type + "/" + diag.name + "/"
        dp.add_string(name_path + "type", diag.type)
        dp.add_string(name_path + "quantity", diag.quantity)
        dp.add_size_t(name_path + "flush_every", diag.flush_every)
        dp.add_array_as_vector(name_path + "write_timestamps", diag.write_timestamps)
        dp.add_array_as_vector(
            name_path + "elapsed_timestamps", diag.elapsed_timestamps
        )
        dp.add_array_as_vector(
            name_path + "compute_timestamps", diag.compute_timestamps
        )

        dp.add_size_t(name_path + "n_attributes", len(diag.attributes))
        for attr_idx, attr_key in enumerate(diag.attributes):
            dp.add_string(name_path + f"attribute_{attr_idx}_key", attr_key)
            if attr_key == "heat_capacity_ratio":
                dp.add_double(
                    name_path + f"attribute_{attr_idx}_value",
                    diag.attributes[attr_key],
                )
            else:
                dp.add_string(
                    name_path + f"attribute_{attr_idx}_value",
                    diag.attributes[attr_key],
                )

    if len(sim.diagnostics) > 0:
        if sim.diag_options is not None and "format" in sim.diag_options:
            dp.add_string(diag_path + "format", sim.diag_options["format"])

        if sim.diag_options is not None and "options" in sim.diag_options:
            dp.add_string(diag_path + "filePath", sim.diag_options["options"]["dir"])
            if "mode" in sim.diag_options["options"]:
                dp.add_string(diag_path + "mode", sim.diag_options["options"]["mode"])
            if "fine_dump_lvl_max" in sim.diag_options["options"]:
                dp.add_int(
                    diag_path + "fine_dump_lvl_max",
                    sim.diag_options["options"]["fine_dump_lvl_max"],
                )
            if "allow_emergency_dumps" in sim.diag_options["options"]:
                dp.add_bool(
                    diag_path + "allow_emergency_dumps",
                    sim.diag_options["options"]["allow_emergency_dumps"],
                )
        else:
            dp.add_string(diag_path + "filePath", "phare_output")


def populate_restarts(dp, sim, serialized_sim):
    if sim.restart_options is None:
        return

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
            raise ValueError(f"PHARE restart file not found for time {restart_time}")

        deserialized_simulation = deserialize_sim(
            _serialized_simulation_string(restart_file_load_path)
        )
        if not sim.is_restartable_compared_to(deserialized_simulation):
            raise ValueError(
                "deserialized Restart simulation is incompatible with configured simulation parameters"
            )

        dp.add_vector_int(
            restarts_path + "restart_ids", _patch_data_ids(restart_file_load_path)
        )
        dp.add_string(restarts_path + "loadPath", restart_file_load_path)
        dp.add_double(restarts_path + "restart_time", restart_time)

    if "mode" in restart_options:
        dp.add_string(restarts_path + "mode", restart_options["mode"])

    dp.add_string(restarts_path + "filePath", restart_file_path)

    if "elapsed_timestamps" in restart_options:
        dp.add_array_as_vector(
            restarts_path + "elapsed_timestamps",
            restart_options["elapsed_timestamps"],
        )

    if "timestamps" in restart_options:
        dp.add_array_as_vector(
            restarts_path + "write_timestamps", restart_options["timestamps"]
        )

    dp.add_string(restarts_path + "serialized_simulation", serialized_sim)


def _populate_refinement_boxes(dp, refinement_boxes):
    dp.add_int(
        "simulation/AMR/refinement/boxes/nbr_levels/", len(refinement_boxes.keys())
    )
    for level, boxes in refinement_boxes.items():
        level_path = "simulation/AMR/refinement/boxes/" + level + "/"
        dp.add_int(level_path + "nbr_boxes/", int(len(boxes)))
        for box_i, box in enumerate(boxes):
            box_id = "B" + str(box_i)
            lower = box.lower
            upper = box.upper
            box_lower_path_x = box_id + "/lower/x/"
            box_upper_path_x = box_id + "/upper/x/"
            dp.add_int(level_path + box_lower_path_x, lower[0])
            dp.add_int(level_path + box_upper_path_x, upper[0])
            if len(lower) >= 2:
                box_lower_path_y = box_id + "/lower/y/"
                box_upper_path_y = box_id + "/upper/y/"
                dp.add_int(level_path + box_lower_path_y, lower[1])
                dp.add_int(level_path + box_upper_path_y, upper[1])
                if len(lower) == 3:
                    box_lower_path_z = box_id + "/lower/z/"
                    box_upper_path_z = box_id + "/upper/z/"
                    dp.add_int(level_path + box_lower_path_z, lower[2])
                    dp.add_int(level_path + box_upper_path_z, upper[2])


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
