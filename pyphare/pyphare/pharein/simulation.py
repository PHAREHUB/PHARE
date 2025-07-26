import os
from copy import deepcopy

import numpy as np

from ..core import box as boxm
from ..core import phare_utilities
from ..core.box import Box
from . import global_vars

# ------------------------------------------------------------------------------


def supported_dimensions():
    return [1, 2]


def compute_dimension(cells):
    return len(cells)


# ------------------------------------------------------------------------------


def de_numpify_recursive(owner, key):
    obj = getattr(owner, key)
    if obj is None:
        return
    if isinstance(obj, np.ndarray):
        object.__setattr__(owner, key, obj.tolist())
    elif hasattr(obj, "__dict__"):
        for k in obj.__dict__:
            de_numpify_recursive(obj, k)


def de_numpify_simulation(sim):
    assert isinstance(sim, Simulation)
    for k in sim.__dict__:
        de_numpify_recursive(sim, k)
    return sim


def re_numpify_recursive(owner, key):
    obj = getattr(owner, key)
    if obj is None:
        return
    if isinstance(obj, list) and len(obj) and not isinstance(obj[0], str):
        object.__setattr__(owner, key, np.array(obj))
    elif hasattr(obj, "__dict__"):
        for k in obj.__dict__:
            re_numpify_recursive(obj, k)


def re_numpify_simulation(sim):
    assert isinstance(sim, Simulation)
    for k in sim.__dict__:
        re_numpify_recursive(sim, k)
    return sim


def check_domain(**kwargs):
    """
    check parameters relative to the domain and raise exceptions if they are not valid
    return 'dl' and 'cells'
    """
    dom_and_dl = "domain_size" in kwargs and "dl" in kwargs and "cells" not in kwargs
    cells_and_dl = "cells" in kwargs and "dl" in kwargs and "domain_size" not in kwargs
    dom_and_cells = "cells" in kwargs and "domain_size" in kwargs and "dl" not in kwargs

    if not (dom_and_dl or dom_and_cells or cells_and_dl):
        raise ValueError(
            "Error: Specify either 'domain_size' and 'dl' or 'cells' and 'dl' or 'domain_size' and 'dl'"
        )

    if dom_and_dl:
        try:
            phare_utilities.check_iterables(kwargs["domain_size"], kwargs["dl"])
            phare_utilities.check_equal_size(kwargs["domain_size"], kwargs["dl"])
        except ValueError:
            raise ValueError(
                "Error: 'domain_size' and 'dl' must have the same dimension"
            )

        domain_size = phare_utilities.listify(kwargs["domain_size"])
        dl = phare_utilities.listify(kwargs["dl"])

        cells = [
            int(dom_size / float(mesh_size))
            for dom_size, mesh_size in zip(domain_size, dl)
        ]

        for s, d in zip(domain_size, dl):
            if s < d:
                raise ValueError("Error: mesh size should be smaller than domain_size")

        dl = [dom_size / n for dom_size, n in zip(domain_size, cells)]

    elif cells_and_dl:
        try:
            phare_utilities.check_iterables(kwargs["cells"], kwargs["dl"])
            phare_utilities.check_equal_size(kwargs["cells"], kwargs["dl"])
        except ValueError:
            raise ValueError(
                "Error: 'cells' and 'dl' must have the same dimension ({} and {})".format(
                    kwargs["cells"], kwargs["dl"]
                )
            )

        dl = phare_utilities.listify(kwargs["dl"])
        cells = phare_utilities.listify(kwargs["cells"])

    elif dom_and_cells:
        try:
            phare_utilities.check_iterables(kwargs["cells"], kwargs["domain_size"])
            phare_utilities.check_equal_size(kwargs["cells"], kwargs["domain_size"])
        except ValueError:
            raise ValueError(
                "Error: 'cells' and 'domain_size' must have the same dimension"
            )

        cells = phare_utilities.listify(kwargs["cells"])
        domain_size = phare_utilities.listify(kwargs["domain_size"])
        dl = [dom_size / float(n) for dom_size, n in zip(domain_size, cells)]

    for n, d in zip(cells, dl):
        if n != 0 and n < 10:
            raise ValueError(
                "Error: number of cells in non-invariant directions should be >= 10"
            )
        if n == 0 and d != 0:
            raise ValueError("Error: dl should be 0 for invariant directions")
        if n < 0:
            raise ValueError("Error: number of cells must be >= 0")
        if d < 0:
            raise ValueError("Error: mesh size must be >= 0")

    return dl, cells


# ------------------------------------------------------------------------------


def check_time(**kwargs):
    """
    check parameters relative to simulation duration and time resolution and raise exception if invalids
     return time_step_nbr and time_step
    """
    final_and_dt = (
        "final_time" in kwargs
        and "time_step" in kwargs
        and "time_step_nbr" not in kwargs
    )
    nsteps_and_dt = (
        "time_step_nbr" in kwargs
        and "time_step" in kwargs
        and "final_time" not in kwargs
    )
    final_and_nsteps = (
        "final_time" in kwargs
        and "time_step_nbr" in kwargs
        and "time_step" not in kwargs
    )

    start_time = kwargs.get("restart_options", {}).get("restart_time", 0)

    if final_and_dt:
        time_step_nbr = int(kwargs["final_time"] / kwargs["time_step"])
        time_step = kwargs["final_time"] / time_step_nbr

    elif final_and_nsteps:
        time_step = kwargs["final_time"] / kwargs["time_step_nbr"]
        time_step_nbr = kwargs["time_step_nbr"]

    elif nsteps_and_dt:
        time_step = kwargs["time_step"]
        time_step_nbr = kwargs["time_step_nbr"]

    else:
        raise ValueError(
            "Error: Specify either 'final_time' and 'time_step' or 'time_step_nbr' and 'time_step'"
            + " or 'final_time' and 'time_step_nbr'"
        )

    return (
        time_step_nbr,
        time_step,
        kwargs.get("final_time", start_time + time_step * time_step_nbr),
    )


# ------------------------------------------------------------------------------


def check_interp_order(**kwargs):
    interp_order = kwargs.get("interp_order", 1)

    if interp_order not in [1, 2, 3]:
        raise ValueError("Error: invalid interpolation order. Should be in [1,2,3]")

    return interp_order


# ------------------------------------------------------------------------------


def check_pusher(**kwargs):
    pusher = kwargs.get("particle_pusher", "modified_boris")
    if pusher not in ["modified_boris"]:
        raise ValueError("Error: invalid pusher ({})".format(pusher))
    return pusher


# ------------------------------------------------------------------------------


def check_layout(**kwargs):
    layout = kwargs.get("layout", "yee")
    if layout not in ("yee"):
        raise ValueError("Error: invalid layout ({})".format(layout))
    return layout


# ------------------------------------------------------------------------------


def check_path(**kwargs):
    path = kwargs.get("path", "." + os.path.sep)
    return path


# ------------------------------------------------------------------------------


def check_boundaries(ndim, **kwargs):
    valid_boundary_types = ("periodic",)
    boundary_types = kwargs.get("boundary_types", ["periodic"] * ndim)
    phare_utilities.check_iterables(boundary_types)

    if phare_utilities.none_iterable(boundary_types):
        bc_length = 1
        if boundary_types not in valid_boundary_types:
            raise ValueError(
                "Error: '{}' is not a valid boundary type".format(boundary_types)
            )
        boundary_types = phare_utilities.listify(boundary_types)
    else:
        bc_length = len(boundary_types)
        for bc in boundary_types:
            if bc not in valid_boundary_types:
                raise ValueError("Error: '{}' is not a valid boundary type".format(bc))

    if bc_length != ndim:
        raise ValueError(
            "Error- boundary_types should have length {} and is of length {}".format(
                ndim, bc_length
            )
        )

    return boundary_types


# ------------------------------------------------------------------------------


# See: https://github.com/PHAREHUB/PHARE/wiki/exactSplitting
# This should match possibleSimulators() in meta_utilities.h
valid_refined_particle_nbr = {
    # dim : {interp : [valid list refined_particle_nbr]}
    1: {1: [2, 3], 2: [2, 3, 4], 3: [2, 3, 4, 5]},
    2: {1: [4, 5, 8, 9], 2: [4, 5, 8, 9, 16], 3: [4, 5, 8, 9, 25]},
    3: {
        1: [6, 7, 8, 9, 12, 13, 14, 15, 18, 19, 20, 21, 26, 27],
        2: [6, 7, 8, 9, 12, 13, 14, 15, 18, 19, 20, 21, 26, 27, 64],
        3: [6, 7, 8, 9, 12, 13, 14, 15, 18, 19, 20, 21, 26, 27, 125],
    },
}  # Default refined_particle_nbr per dim/interp is considered index 0 of list


def check_refined_particle_nbr(ndim, **kwargs):
    interp = kwargs["interp_order"]
    refined_particle_nbr = kwargs.get(
        "refined_particle_nbr", valid_refined_particle_nbr[ndim][interp][0]
    )

    if refined_particle_nbr not in valid_refined_particle_nbr[ndim][interp]:
        raise ValueError(
            "Invalid split particle number, valid values for dim({}) ".format(ndim)
            + "interp({}) include {}".format(
                interp, valid_refined_particle_nbr[ndim][interp]
            )
        )

    return refined_particle_nbr


# ------------------------------------------------------------------------------


def check_origin(ndim, **kwargs):
    origin = kwargs.get("origin", [0.0] * ndim)
    return origin


# ------------------------------------------------------------------------------


def as_list_per_level(refinement_boxes):
    """
    accepts various formats of boxes.
    returns {"L0" : [Box(),Box()]}
    """

    list_per_level = {}

    for level_key, boxes in refinement_boxes.items():
        if isinstance(level_key, str):
            assert level_key[0] == "L"
        else:
            level_key = "L" + str(level_key)

        if isinstance(boxes, list) and all([isinstance(box, Box) for box in boxes]):
            list_per_level[level_key] = boxes

        elif isinstance(boxes, dict):
            list_per_level[level_key] = []

            if all([isinstance(val, Box) for key, val in boxes.items()]):
                for box_id, box in boxes.items():
                    list_per_level[level_key] += [box]

            elif all([isinstance(val, (tuple, list)) for key, val in boxes.items()]):
                for box_id, box_coords in boxes.items():
                    box_coords = [list(box_coord) for box_coord in box_coords]
                    list_per_level[level_key] += [Box(box_coords[0], box_coords[1])]

    return list_per_level


def check_refinement_boxes(ndim, **kwargs):
    """
    returns tuple ( { "L0" : [Boxes]}, max_nbr_levels)
    """

    refinement_boxes = kwargs.get("refinement_boxes", None)

    if (
        refinement_boxes is None or len(refinement_boxes) == 0
    ):  # SAMRAI errors if empty dict
        return None, 1  # default max level number to 1 if no boxes

    smallest_patch_size = kwargs.get("smallest_patch_size")
    refinement_ratio = kwargs["refinement_ratio"]
    nesting_buffer = kwargs["nesting_buffer"]

    refinement_boxes = as_list_per_level(refinement_boxes)
    domain_box = Box([0] * ndim, np.asarray(kwargs["cells"]) - 1)

    boxes_per_level = {0: [domain_box]}

    for level_key, boxes in refinement_boxes.items():
        ilvl = int(level_key[1:])
        refined_level_number = ilvl + 1
        boxes_per_level[refined_level_number] = [
            boxm.refine(box, refinement_ratio) for box in boxes
        ]

        if len(boxes) == 0:
            raise ValueError("Error - missing refinement boxes")

        for box_idx, ref_box in enumerate(boxes):
            for cmp_box in boxes[box_idx + 1 :]:
                if ref_box * cmp_box is not None:
                    raise ValueError(
                        f"Error: Box({ref_box}) overlaps with Box({cmp_box})"
                    )

        for box in boxes:
            refined_coarser_boxes = boxes_per_level[ilvl]

            coarse_nCell = 0
            for refined_coarser in refined_coarser_boxes:
                intersection = box * boxm.shrink(refined_coarser, nesting_buffer)
                coarse_nCell += 0 if intersection is None else intersection.nCells()

            if coarse_nCell != box.nCells():
                raise ValueError(
                    f"Box({box}) is incompatible with coarser boxes({refined_coarser_boxes}) and nest_buffer({nesting_buffer})"
                )

            if box.ndim != ndim:
                raise ValueError(f"Box({box}) has incorrect dimensions for simulation")
            for length in boxm.refine(box, refinement_ratio).shape:
                if (length < smallest_patch_size).any():
                    raise ValueError(
                        "Invalid box incompatible with smallest_patch_size"
                    )

    return refinement_boxes, len(refinement_boxes.items()) + 1


# ------------------------------------------------------------------------------


def check_patch_size(ndim, **kwargs):
    def get_max_ghosts():
        from ..core.gridlayout import GridLayout

        grid = GridLayout()
        return max(
            grid.nbrGhosts(kwargs["interp_order"], x) for x in ["primal", "dual"]
        )

    max_ghosts = get_max_ghosts()
    small_invalid_patch_size = phare_utilities.np_array_ify(max_ghosts, ndim)
    largest_patch_size = kwargs.get("largest_patch_size", None)

    # to prevent primal ghost overlaps of non adjacent patches, we need smallest_patch_size+=1
    smallest_patch_size = phare_utilities.np_array_ify(max_ghosts, ndim) + 1
    if "smallest_patch_size" in kwargs and kwargs["smallest_patch_size"] is not None:
        smallest_patch_size = phare_utilities.np_array_ify(
            kwargs["smallest_patch_size"], ndim
        )

    cells = phare_utilities.np_array_ify(kwargs["cells"])

    if smallest_patch_size.size != ndim:
        raise ValueError(
            f"Error: smallest_patch_size({smallest_patch_size.size}) must be size {ndim}"
        )

    if (smallest_patch_size <= small_invalid_patch_size).any():
        raise ValueError(
            "Error - smallest_patch_size cannot be <= " + str(small_invalid_patch_size)
        )

    if (smallest_patch_size > cells).any():
        raise ValueError(
            "Error - smallest_patch_size should be less than nbr of cells in all directions"
        )

    if largest_patch_size is not None:
        largest_patch_size = phare_utilities.np_array_ify(largest_patch_size, ndim)

        if largest_patch_size.size != ndim:
            raise ValueError(
                f"Error: largest_patch_size({largest_patch_size.size}) must be size {ndim}"
            )

        if (largest_patch_size > cells).any():
            raise ValueError(
                "Error - largest_patch_size should be less than nbr of cells in all directions"
            )

        if (largest_patch_size <= 0).any():
            raise ValueError("Error - largest_patch_size cannot be <= 0")

        if (largest_patch_size < smallest_patch_size).any():
            raise ValueError(
                "Error - largest_patch_size and smallest_patch_size are incompatible"
            )

    return largest_patch_size, smallest_patch_size


# ------------------------------------------------------------------------------


def check_directory(directory, key):
    directory = directory.rstrip(os.path.sep)  # trim trailing slashes

    if os.path.exists(directory):
        if os.path.isfile(directory):
            raise ValueError(f"Error: Simulation {key} dir exists as a file.")

        if not os.access(directory, os.R_OK | os.W_OK | os.X_OK):
            raise ValueError(
                f"Directory ({directory}) for {key} does not have the correct permissions"
            )

    try:
        os.makedirs(directory, exist_ok=True)
        if not os.path.exists(directory):
            raise ValueError("1. Creation of the directory %s failed" % directory)
    except FileExistsError:
        raise ValueError("Creation of the directory %s failed" % directory)

    return directory


# diag_options = {"format":"phareh5", "options": {"dir": "phare_ouputs/"}}
def check_diag_options(**kwargs):
    diag_options = kwargs.get("diag_options", None)
    formats = ["phareh5"]
    if diag_options is not None and "format" in diag_options:
        if diag_options["format"] not in formats:
            raise ValueError("Error - diag_options format is invalid")
        if "options" in diag_options and "dir" in diag_options["options"]:
            diag_options["options"]["dir"] = check_directory(
                diag_options["options"]["dir"], "diagnostics"
            )
        valid_modes = ["overwrite"]
        if "options" in diag_options and "mode" in diag_options["options"]:
            mode = diag_options["options"]["mode"]
            if mode not in valid_modes:
                raise ValueError(
                    f"Invalid diagnostics mode {mode}, valid modes are {valid_modes}"
                )
    return diag_options


def check_restart_options(**kwargs):
    valid_keys = ["dir", "elapsed_timestamps", "timestamps", "mode", "restart_time"]
    restart_options = kwargs.get("restart_options", None)

    if restart_options is not None:
        for key in restart_options.keys():
            if key not in valid_keys:
                raise ValueError(
                    f"invalid option ({key}), valid options are {valid_keys}"
                )

        valid_modes = ["conserve", "overwrite"]
        if "mode" not in restart_options:
            raise ValueError(f"Restart mode not set, valid modes are {valid_modes}")

        if "dir" in restart_options:
            restart_options["dir"] = check_directory(restart_options["dir"], "restarts")

        mode = restart_options["mode"]
        if mode not in valid_modes:
            raise ValueError(
                f"Invalid restart mode {mode}, valid modes are {valid_modes}"
            )

    return restart_options


def validate_restart_options(sim):
    import pyphare.pharein.restarts as restarts

    if sim.restart_options is not None:
        restarts.validate(sim)


def check_refinement(**kwargs):
    return kwargs.get("refinement", "boxes")


def check_nesting_buffer(ndim, **kwargs):
    nesting_buffer = phare_utilities.np_array_ify(kwargs.get("nesting_buffer", 0), ndim)

    if nesting_buffer.size != ndim:
        raise ValueError(f"Error: nesting_buffer must be size {ndim}")

    if (nesting_buffer < 0).any():
        raise ValueError(f"Error: nesting_buffer({nesting_buffer}) cannot be negative")

    smallest_patch_size = kwargs["smallest_patch_size"]
    largest_patch_size = kwargs["largest_patch_size"]

    if (nesting_buffer > smallest_patch_size / 2).any():
        raise ValueError(
            f"Error: nesting_buffer({nesting_buffer})"
            + f"cannot be larger than half the smallest_patch_size ({smallest_patch_size})"
        )

    cells = np.asarray(kwargs["cells"])

    if (
        largest_patch_size is not None
        and (nesting_buffer > (cells - largest_patch_size)).any()
    ):
        raise ValueError(
            f"Error: nesting_buffer({nesting_buffer})"
            + f"incompatible with number of domain cells ({cells}) and largest_patch_size({largest_patch_size})"
        )

    elif (nesting_buffer > cells).any():
        raise ValueError(
            f"Error: nesting_buffer({nesting_buffer}) incompatible with number of domain cells ({cells})"
        )

    return nesting_buffer


def check_optional_keywords(**kwargs):
    extra = []

    if check_refinement(**kwargs) != "boxes":
        extra += ["max_nbr_levels"]

    return extra


def check_resistivity(**kwargs):
    resistivity = kwargs.get("resistivity", 0.0)
    if resistivity < 0.0:
        raise ValueError("Error: resistivity should not be negative")

    return resistivity


def check_hyper_resistivity(**kwargs):
    hyper_resistivity = kwargs.get("hyper_resistivity", 0.0001)
    if hyper_resistivity < 0.0:
        raise ValueError("Error: hyper_resistivity should not be negative")

    hyper_mode = kwargs.get("hyper_mode", "constant")
    if hyper_mode not in ["constant", "spatial"]:
        raise ValueError("Error: hyper_mode should be 'constant' or 'spatial'")

    return hyper_resistivity, hyper_mode


def check_clustering(**kwargs):
    valid_keys = ["berger", "tile"]
    clustering = kwargs.get("clustering", "berger")
    if clustering not in valid_keys:
        raise ValueError(
            f"Error: clustering type is not supported, supported types are {valid_keys}"
        )
    return clustering


def check_max_mhd_level(**kwargs):
    max_mhd_level = kwargs.get("max_mhd_level", 0)

    if max_mhd_level > kwargs["max_nbr_levels"]:
        raise ValueError(
            f"Error: max_mhd_level({max_mhd_level}) should be less or equal to max_nbr_levels({kwargs['max_nbr_levels']})"
        )

    return max_mhd_level

def check_model_options(**kwargs):
    model_options = kwargs.get("model_options", None)

    if model_options is None:
        return None

    valid_options = {"MHDModel", "HybridModel"}

    if not set(model_options).issubset(valid_options):
        raise ValueError(f"Invalid model options: {model_options}. Allowed values are {valid_options}.")

    return model_options

def check_mhd_constants(**kwargs):
    gamma = kwargs.get("gamma", 5.0 / 3.0)
    eta = kwargs.get("eta", 0.0)
    nu = kwargs.get("nu", 0.0)

    return gamma, eta, nu

def check_mhd_terms(**kwargs):
    hall = kwargs.get("hall", False)
    res = kwargs.get("res", False)
    hyper_res = kwargs.get("hyper_res", False)

    return hall, res, hyper_res

def check_mhd_parameters(**kwargs):
    reconstruction = kwargs.get("reconstruction", "")
    limiter = kwargs.get("limiter", "")
    riemann = kwargs.get("riemann", "")
    mhd_timestepper = kwargs.get("mhd_timestepper", "")

    return reconstruction, limiter, riemann, mhd_timestepper


# ------------------------------------------------------------------------------


def checker(func):
    def wrapper(simulation_object, **kwargs):
        accepted_keywords = [
            "domain_size",
            "cells",
            "dl",
            "particle_pusher",
            "final_time",
            "time_step",
            "time_step_nbr",
            "layout",
            "interp_order",
            "origin",
            "boundary_types",
            "refined_particle_nbr",
            "path",
            "nesting_buffer",
            "diag_export_format",
            "refinement_boxes",
            "refinement",
            "tagging_threshold",
            "clustering",
            "smallest_patch_size",
            "largest_patch_size",
            "diag_options",
            "resistivity",
            "hyper_resistivity",
            "hyper_mode",
            "strict",
            "restart_options",
            "tag_buffer",
            "description",
            "dry_run",
            "write_reports",

            "max_mhd_level",
            "model_options",
            "gamma",
            "eta",
            "nu",
            "hall",
            "res",
            "hyper_res",
            "reconstruction",
            "limiter",
            "riemann",
            "mhd_timestepper",
        ]

        accepted_keywords += check_optional_keywords(**kwargs)

        wrong_kwds = phare_utilities.not_in_keywords_list(accepted_keywords, **kwargs)
        if len(wrong_kwds) > 0:
            raise ValueError("Error: invalid arguments - " + " ".join(wrong_kwds))

        dl, cells = check_domain(**kwargs)

        kwargs["strict"] = kwargs.get("strict", False)

        kwargs["dl"] = dl
        kwargs["cells"] = cells
        kwargs["refinement_ratio"] = 2

        kwargs["description"] = kwargs.get("description", None)

        kwargs["clustering"] = check_clustering(**kwargs)

        time_step_nbr, time_step, final_time = check_time(**kwargs)
        kwargs["time_step_nbr"] = time_step_nbr
        kwargs["time_step"] = time_step
        kwargs["final_time"] = final_time

        kwargs["interp_order"] = check_interp_order(**kwargs)
        kwargs["refinement_ratio"] = 2

        kwargs["particle_pusher"] = check_pusher(**kwargs)
        kwargs["layout"] = check_layout(**kwargs)
        kwargs["path"] = check_path(**kwargs)

        ndim = compute_dimension(cells)
        kwargs["diag_options"] = check_diag_options(**kwargs)
        kwargs["restart_options"] = check_restart_options(**kwargs)

        kwargs["boundary_types"] = check_boundaries(ndim, **kwargs)
        kwargs["origin"] = check_origin(ndim, **kwargs)

        kwargs["refined_particle_nbr"] = check_refined_particle_nbr(ndim, **kwargs)
        kwargs["diag_export_format"] = kwargs.get("diag_export_format", "hdf5")
        assert kwargs["diag_export_format"] in ["hdf5"]  # only hdf5 supported for now

        largest, smallest = check_patch_size(ndim, **kwargs)
        kwargs["smallest_patch_size"] = smallest
        kwargs["largest_patch_size"] = largest

        kwargs["nesting_buffer"] = check_nesting_buffer(ndim, **kwargs)

        kwargs["tag_buffer"] = kwargs.get("tag_buffer", 1)

        kwargs["refinement"] = check_refinement(**kwargs)
        if kwargs["refinement"] == "boxes":
            (
                kwargs["refinement_boxes"],
                kwargs["max_nbr_levels"],
            ) = check_refinement_boxes(ndim, **kwargs)
        else:
            kwargs["max_nbr_levels"] = kwargs.get("max_nbr_levels", None)
            assert kwargs["max_nbr_levels"] is not None  # this needs setting otherwise
            kwargs["refinement_boxes"] = None
            kwargs["tagging_threshold"] = kwargs.get("tagging_threshold", 0.1)

        kwargs["resistivity"] = check_resistivity(**kwargs)

        nu, hyper_mode = check_hyper_resistivity(**kwargs)
        kwargs["hyper_resistivity"] = nu
        kwargs["hyper_mode"] = hyper_mode

        kwargs["dry_run"] = kwargs.get(
            "dry_run", os.environ.get("PHARE_DRY_RUN", "0") == "1"
        )
        kwargs["write_reports"] = kwargs.get(  # on by default except for tests
            "write_reports", os.environ.get("PHARE_TESTING", "0") != "1"
        )

        kwargs["max_mhd_level"] = check_max_mhd_level(**kwargs)

        kwargs["model_options"] = check_model_options(**kwargs)

        gamma, eta, nu = check_mhd_constants(**kwargs)
        kwargs["gamma"] = gamma
        kwargs["eta"] = eta
        kwargs["nu"] = nu

        hall, res, hyper_res = check_mhd_terms(**kwargs)
        kwargs["hall"] = hall
        kwargs["res"] = res
        kwargs["hyper_res"] = hyper_res

        reconstruction, limiter, riemann, mhd_timestepper = check_mhd_parameters(**kwargs)
        kwargs["reconstruction"] = reconstruction
        kwargs["limiter"] = limiter
        kwargs["riemann"] = riemann
        kwargs["mhd_timestepper"] = mhd_timestepper

        return func(simulation_object, **kwargs)

    return wrapper


# ------------------------------------------------------------------------------


class Simulation(object):
    """

    **Usage example:**

    This declares a 2D simulation of 100x100 cells
    with an isotropic mesh size of 0.2. The simulation will run
    1000 time steps of dt=0.001 and so will stop at t=1

    The simulation will run with adaptive mesh refinement and evolve up to 3 levels,
    i.e. a base coarsest mesh with up to 2 additional refinement levels.

    Diagnostics, if declared, will be saved as native PHARE HDF5 (default format)
    in the directory `diag_outputs` and will overwrite any existing h5 files there.

    The resistivity and hyper-resistivity are set to constant custom values.

        .. code-block:: python

            from pyphare.pharein import Simulation



            Simulation(
                       # time evolution parameters
                       time_step_nbr=1000,
                       time_step=0.001,

                       # Domain parameters
                       cells=(100,100),
                       dl=(0.2, 0.2),

                       # AMR parameters
                       max_nbr_levels=3,

                        # general physics parameters
                       hyper_resistivity=0.001,
                       resistivity=0.001,

                       # diagnostics and restart parameters
                       diag_options={"format": "phareh5",
                                     "options": {"dir": diag_outputs,
                                                 "mode":"overwrite"}},
                      )


    **Time parameters:**

        * **final_time** (``float``), final simulation time. Use with time_step OR time_step_nbr
        * **time_step** (``float``), simulation time step. Use with time_step_nbr OR final_time
        * **time_step_nbr** (``int``), number of time step to perform. Use with final_time OR time_step


        .. code-block:: python

            Simulation(
                       time_step_nbr=1000,
                       time_step=0.001,

                    # ....
                    )

        is equivalent to:

        .. code-block:: python

            Simulation(
                       time_step_nbr=1000,
                       final_time=1.,

                    # ....
                    )

        or even to:

        .. code-block:: python

            Simulation(
                       time_step=0.001,
                       final_time=1.,

                    # ....
                    )


    **Domain/grid parameters:**

        The number of dimensions of the simulation is deduced from the length
        of these parameters.

        * **dl** (``float``, ``tuple``) grid spacing dx, (dx,dy) or (dx,dy,dz) in 1D, 2D or 3D.
          A single float value in 2D or 3D are assumed equal for each dimension.
          Use this parameter with either `domain_size` OR `cells`

        * **domain_size** (``float`` or ``tuple``)
          size of the physical domain Lx, (Lx,Ly), (Lx,Ly,Lz) in 1D, 2D or 3D
          Single float is assumed equal for each direction. Use this parameter with either `cells` or `dl`

        * **cells** (``int`` or ``tuple``)
          number of cells nx or (nx, ny) or (nx, ny, nz) in 1, 2 and 3D.
          Single int is assumed equal for each direction
          Use this parameter with either `domain_size` or `dl`

        **Expert parameters:**

        These parameters are more advanced, modify them at your own risk

            * **layout** (``str``), layout of the physical quantities on the mesh (default = "yee")
            * **origin** (``int`` or ``tuple``), origin of the physical domain, (default (0,0,0) in 3D)


    For instance:

        .. code-block:: python

            Simulation(
                       dl=(0.1, 0.1),
                       domain_size=(100,100),

                    # ....
                    )

        is equivalent to:


        .. code-block:: python

            Simulation(
                       cells=(1000, 1000),
                       domain_size=(100,100),

                    # ....
                    )

        or even to:

        .. code-block:: python

            Simulation(
                       cells=(1000, 1000),
                       dl=(0.1,0.1),

                    # ....
                    )

    **Macro-particle parameters:**

        * **interp_order** (``int``), 1, 2 or 3 (default=1) particle b-spline order
        * **particle_pusher** (``str``), algo to push particles (default = "modifiedBoris")


    **Diagnostics output parameters:**

        * **diag_options** (``dict``)
            * **path** (``str``) path for outputs (default : './')
            * **diag_export_format** (``str``) format of the output diagnostics (default= "phareh5")
            * **mode** (``str``) mode of the output diagnostics (default= "overwrite" will write over existing files)



    **Adaptive Mesh Refinement (AMR) parameters:**

        * **max_nbr_levels** (``int``), default=1, max number of levels in the hierarchy. Used if no `refinement_boxes` are set
        * **tag_buffer** (``int``), default=1, value representing the number of cells by which tagged cells are buffered before clustering into boxes. The larger `tag_buffer`, the wider refined regions will be around tagged cells.
        * **clustering** (``str``), {"berger" (default), "tile"}, type of clustering to use for AMR. `tile` results in wider patches, less artifacts and better scalability

        **Expert parameters:**

        These parameters are more advanced, modify them at your own risk

        * **refined_particle_nbr** (``ìnt``), number of refined particle per coarse particle.
        * **nesting_buffer** (``ìnt``), default=0 minimum gap in coarse cells from the border of a level and any refined patch border
        * **refinement_boxes**, default=None, {"L0":{"B0":[(lox,loy,loz),(upx,upy,upz)],...,"Bi":[(),()]},..."Li":{B0:[(),()]}}
        * **smallest_patch_size** (``int`` or ``tuple``), minimum number of cells in a patch in each direction
          This parameter cannot be smaller than the number of field ghost nodes
        * **largest_patch_size** (``int`` or ``tuple``), maximum size of a patch in each direction



    **Restart parameters:**

    These parameters are used to restart a simulation from a previous state and dump chekpoints.

        * **restart_options** (``dict``)
            * **dir** (``str``) path for restart files (default : './')
            * **mode** (``str``) mode of the restart files
                * "conserve"  - (default), will conserve existing files
                * "overwrite" - will overwrite existing files

            * **restart_time** (``float``) time at which to restart the simulation (default=0)
            * **timestamps** (``list``) list of timestamps at which to restart the simulation

    Misc:

        * **description** (``string``), [default=None] arbitrary string for per simulation context - injected in output files when feasible
        * **strict** (``bool``), turns warnings into errors (default False)
        * **resistivity** (``float``), resistivity value (default=0.0)
        * **hyper-resistivity** (``float``), hyper-resistivity value (default=0.0)
        * **boundary_types** (``str`` or ``tuple``) type of boundary conditions (default is "periodic" for each direction)

    """

    @checker
    def __init__(self, **kwargs):
        if global_vars.sim is not None:
            raise RuntimeError("simulation is already created")

        global_vars.sim = self

        for k, v in kwargs.items():
            object.__setattr__(self, k, v)

        self.ndim = compute_dimension(self.cells)

        self.diagnostics = {}
        self.uniform_model = None
        self.maxwellian_fluid_model = None
        self.mhd_model = None
        self.electrons = None
        self.load_balancer = None

        # hard coded in C++ MultiPhysicsIntegrator::getMaxFinerLevelDt
        self.nSubcycles = 4
        self.stepDiff = 1 / self.nSubcycles

        levelNumbers = list(range(self.max_nbr_levels))
        self.level_time_steps = [
            self.time_step * (self.stepDiff ** (ilvl)) for ilvl in levelNumbers
        ]
        self.level_step_nbr = [
            self.nSubcycles ** levelNumbers[ilvl] * self.time_step_nbr
            for ilvl in levelNumbers
        ]
        validate_restart_options(self)

    def simulation_domain(self):
        return [dl * n + ori for dl, n, ori in zip(self.dl, self.cells, self.origin)]

    def within_simulation_duration(self, time_period):
        return time_period[0] >= 0 and time_period[1] < self.time_step_nbr

    def within_simulation_extent(self, extent):
        domain = self.simulation_domain()
        if len(extent) == 2:
            # 1D case
            return extent[0] >= domain[0] and extent[1] <= domain[1]
        raise NotImplementedError("Error: 2D and 3D not implemented yet")

    def restart_file_path(self):
        assert self.restart_options is not None
        return self.restart_options.get("dir", "phare_outputs")

    def start_time(self):
        if self.is_from_restart():
            return self.restart_options["restart_time"]
        return 0

    def is_from_restart(self):
        return (
            self.restart_options is not None and "restart_time" in self.restart_options
        )

    def __getattr__(
        self, name
    ):  # stops pylint complaining about dynamic class attributes not existing
        ...

    # dill serialization uses __getattr__
    #  but without these it causes errors
    def __getstate__(self):
        return vars(self)

    def __setstate__(self, state):
        vars(self).update(state)

    # ------------------------------------------------------------------------------

    def add_diagnostics(self, diag):
        """

        :meta private:

        """
        if diag.name in self.diagnostics:
            raise ValueError(
                "Error: diagnostics {} already registered".format(diag.name)
            )

        # check whether the spatial extent of the diagnostics is valid, given the domain size
        if diag.extent() is not None and not self.within_simulation_extent(
            diag.extent()
        ):
            raise RuntimeError("Error: invalid diagnostics spatial extent")

        self.diagnostics[diag.name] = diag

    # ------------------------------------------------------------------------------

    def is_restartable_compared_to(self, sim):
        """

        :meta private:
        """
        # to be considered https://github.com/PHAREHUB/PHARE/issues/666
        check = ["cells", "dl"]
        are_comparable = all([getattr(self, k) == getattr(sim, k)] for k in check)
        return are_comparable

    # ------------------------------------------------------------------------------

    def count_diagnostics(self, type_name):
        """

        :meta private:
        """
        return len(
            [
                diag
                for diagname, diag in self.diagnostics.items()
                if diag.type == type_name
            ]
        )

    # ------------------------------------------------------------------------------

    def set_uniform_model(self, mhd_model):
        """

        :meta private:
        """
        self.uniform_model = mhd_model

    def set_maxwellian_fluid_model(self, maxwellian_fluid_model):
        """

        :meta private:
        """
        self.maxwellian_fluid_model = maxwellian_fluid_model

    def set_mhd_model(self, mhd_model):
        """

        :meta private:
        """
        self.mhd_model = mhd_model

    def set_electrons(self, electrons):
        """

        :meta private:
        """
        self.electrons = electrons


# ------------------------------------------------------------------------------


def serialize(sim):
    """

    :meta private:
    """
    # pickle cannot handle simulation objects
    import codecs

    import dill

    return codecs.encode(dill.dumps(de_numpify_simulation(deepcopy(sim))), "hex")


def deserialize(hex):
    """:meta private:"""
    import codecs

    import dill

    return re_numpify_simulation(dill.loads(codecs.decode(hex, "hex")))
