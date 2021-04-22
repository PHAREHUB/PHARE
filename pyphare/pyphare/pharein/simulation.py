import os
import numpy as np

from ..core import phare_utilities
from . import global_vars
from ..core import box as boxm
from ..core.box import Box


def supported_dimensions():
    return [1]


def compute_dimension(cells):
    return len(cells)


# ------------------------------------------------------------------------------


def check_domain(**kwargs):
    """
    check parameters relative to the domain and raise exceptions if they are not valid
    return 'dl' and 'cells'
    """
    dom_and_dl = ('domain_size' in kwargs and 'dl' in kwargs and 'cells' not in kwargs)
    cells_and_dl = ('cells' in kwargs and 'dl' in kwargs and 'domain_size' not in kwargs)
    dom_and_cells = ('cells' in kwargs and 'domain_size' in kwargs and 'dl' not in kwargs)

    if not (dom_and_dl or dom_and_cells or cells_and_dl):
        raise ValueError("Error: Specify either 'domain_size' and 'dl' or 'cells' and 'dl' or 'domain_size' and 'dl'")

    if dom_and_dl:
        try:
            phare_utilities.check_iterables(kwargs['domain_size'], kwargs['dl'])
            phare_utilities.check_equal_size(kwargs['domain_size'], kwargs['dl'])
        except ValueError:
            raise ValueError("Error: 'domain_size' and 'dl' must have the same dimension")

        domain_size = phare_utilities.listify(kwargs['domain_size'])
        dl = phare_utilities.listify(kwargs['dl'])

        cells = [int(dom_size / float(mesh_size)) for dom_size, mesh_size in zip(domain_size, dl)]

        for s, d in zip(domain_size, dl):
            if s < d:
                raise ValueError("Error: mesh size should be smaller than domain_size")

        dl = [dom_size / n for dom_size, n in zip(domain_size, cells)]

    elif cells_and_dl:
        try:
            phare_utilities.check_iterables(kwargs['cells'], kwargs['dl'])
            phare_utilities.check_equal_size(kwargs['cells'], kwargs['dl'])
        except ValueError:
            raise ValueError("Error: 'cells' and 'dl' must have the same dimension ({} and {})".format(kwargs['cells'], kwargs["dl"]))

        dl = phare_utilities.listify(kwargs['dl'])
        cells = phare_utilities.listify(kwargs['cells'])

    elif dom_and_cells:
        try:
            phare_utilities.check_iterables(kwargs['cells'], kwargs['domain_size'])
            phare_utilities.check_equal_size(kwargs['cells'], kwargs['domain_size'])
        except ValueError:
            raise ValueError("Error: 'cells' and 'domain_size' must have the same dimension")

        cells = phare_utilities.listify(kwargs['cells'])
        domain_size = phare_utilities.listify(kwargs['domain_size'])
        dl = [dom_size / float(n) for dom_size, n in zip(domain_size, cells)]

    for n, d in zip(cells, dl):
        if n != 0 and n < 10:
            raise ValueError("Error: number of cells in non-invariant directions should be >= 10")
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
    final_and_dt = ('final_time' in kwargs and 'time_step' in kwargs and 'time_step_nbr' not in kwargs)
    nsteps_and_dt = ('time_step_nbr' in kwargs and 'time_step' in kwargs and 'final_time' not in kwargs)
    final_and_nsteps = ('final_time' in kwargs and 'time_step_nbr' in kwargs and 'time_step' not in kwargs)

    if final_and_dt:
        time_step_nbr = int(kwargs['final_time'] / kwargs['time_step'])
        time_step = kwargs['final_time'] / time_step_nbr

    elif final_and_nsteps:
        time_step = kwargs['final_time'] / kwargs['time_step_nbr']
        time_step_nbr = kwargs['time_step_nbr']

    elif nsteps_and_dt:
        time_step = kwargs['time_step']
        time_step_nbr = kwargs['time_step_nbr']

    else:
        raise ValueError("Error: Specify either 'final_time' and 'time_step' or 'time_step_nbr' and 'time_step'" + \
                         " or 'final_time' and 'time_step_nbr'")

    return time_step_nbr, time_step, kwargs.get('final_time', time_step * time_step_nbr)


# ------------------------------------------------------------------------------


def check_interp_order(**kwargs):
    interp_order = kwargs.get('interp_order', 1)

    if interp_order not in [1, 2, 3]:
        raise ValueError("Error: invalid interpolation order. Should be in [1,2,3]")

    return interp_order

# ------------------------------------------------------------------------------


def check_pusher(**kwargs):
    pusher = kwargs.get('particle_pusher', 'modified_boris')
    if pusher not in ['modified_boris', 'kirov']:
        raise ValueError('Error: invalid pusher ({})'.format(pusher))
    return pusher


# ------------------------------------------------------------------------------


def check_layout(**kwargs):
    layout = kwargs.get('layout', 'yee')
    if layout not in ('yee'):
        raise ValueError('Error: invalid layout ({})'.format(layout))
    return layout


# ------------------------------------------------------------------------------


def check_path(**kwargs):
    path = kwargs.get('path', '.' + os.path.sep)
    return path

# ------------------------------------------------------------------------------


def check_boundaries(ndim, **kwargs):
    valid_boundary_types = ("periodic",)
    boundary_types = kwargs.get('boundary_types', ['periodic'] * ndim)
    phare_utilities.check_iterables(boundary_types)

    if phare_utilities.none_iterable(boundary_types):
        bc_length = 1
        if boundary_types not in valid_boundary_types:
            raise ValueError("Error: '{}' is not a valid boundary type".format(boundary_types))
        boundary_types = phare_utilities.listify(boundary_types)
    else:
        bc_length = len(boundary_types)
        for bc in boundary_types:
            if bc not in valid_boundary_types:
                raise ValueError("Error: '{}' is not a valid boundary type".format(bc))

    if bc_length != ndim:
        raise ValueError("Error- boundary_types should have length {} and is of length {}".format(ndim, bc_length))

    return boundary_types


# ------------------------------------------------------------------------------


# See: https://github.com/PHAREHUB/PHARE/wiki/exactSplitting
# This should match possibleSimulators() in meta_utilities.h
valid_refined_particle_nbr = {
  # dim : {interp : [valid list refined_particle_nbr]}
  1: {
    1: [2, 3],
    2: [2, 3, 4],
    3: [2, 3, 4, 5]
  },
  2: {
    1: [4, 5, 8, 9],
    2: [4, 5, 8, 9, 16],
    3: [4, 5, 8, 9, 25]
  },
  3: {
    1: [6, 7, 8, 9, 12, 13, 14, 15, 18, 19, 20, 21, 26, 27],
    2: [6, 7, 8, 9, 12, 13, 14, 15, 18, 19, 20, 21, 26, 27, 64],
    3: [6, 7, 8, 9, 12, 13, 14, 15, 18, 19, 20, 21, 26, 27, 125]
  }
} # Default refined_particle_nbr per dim/interp is considered index 0 of list
def check_refined_particle_nbr(ndim, **kwargs):

    interp = kwargs["interp_order"]
    refined_particle_nbr = kwargs.get("refined_particle_nbr", valid_refined_particle_nbr[ndim][interp][0])

    if refined_particle_nbr not in valid_refined_particle_nbr[ndim][interp]:
        raise ValueError("Invalid split particle number, valid values for dim({}) ".format(ndim)
            + "interp({}) include {}".format(interp, valid_refined_particle_nbr[ndim][interp]))

    return refined_particle_nbr



# ------------------------------------------------------------------------------


def check_origin(ndim, **kwargs):
    origin = kwargs.get("origin", [0.] * ndim)
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
            level_key = "L"+str(level_key)

        if isinstance(boxes, list) and all([isinstance(box, Box) for box in boxes]):
          list_per_level[level_key] = boxes

        elif isinstance(boxes, dict):
            list_per_level[level_key] = []

            if all([isinstance(val, Box) for key,val in boxes.items()]):
                for box_id, box in boxes.items():
                    list_per_level[level_key] += [box]

            elif all([isinstance(val, (tuple,list)) for key,val in boxes.items()]):
                for box_id, box_coords in boxes.items():
                    box_coords = [list(box_coord) for box_coord in box_coords]
                    list_per_level[level_key] += [Box(box_coords[0], box_coords[1])]

    return list_per_level

def check_refinement_boxes(ndim, **kwargs):
    """
      returns tuple ( { "L0" : [Boxes]}, max_nbr_levels)
    """

    refinement_boxes = kwargs.get("refinement_boxes", None)

    if refinement_boxes is None or len(refinement_boxes) == 0: # SAMRAI errors if empty dict
        return None, 1 # default max level number to 1 if no boxes

    smallest_patch_size = kwargs.get("smallest_patch_size")
    refinement_ratio = kwargs["refinement_ratio"]
    nesting_buffer = kwargs["nesting_buffer"]

    refinement_boxes = as_list_per_level(refinement_boxes)
    domain_box = Box([0] * ndim, np.asarray(kwargs["cells"]) - 1)

    boxes_per_level = {0: [domain_box]}

    for level_key, boxes in refinement_boxes.items():

        ilvl = int(level_key[1:])
        refined_level_number = ilvl+1
        boxes_per_level[refined_level_number] = [boxm.refine(box, refinement_ratio) for box in boxes]

        if len(boxes) == 0:
            raise ValueError("Error - missing refinement boxes")

        for box_idx, ref_box in enumerate(boxes):
            for cmp_box in boxes[box_idx + 1:]:
                if ref_box * cmp_box != None:
                    raise ValueError(f"Error: Box({ref_box}) overlaps with Box({cmp_box})")

        for box in boxes:
            refined_coarser_boxes = boxes_per_level[ilvl]

            coarse_nCell = 0
            for refined_coarser in refined_coarser_boxes:
                intersection = box * boxm.shrink(refined_coarser, [nesting_buffer] * ndim)
                coarse_nCell += 0 if intersection is None else intersection.nCells()

            if coarse_nCell != box.nCells():
                raise ValueError(f"Box({box}) is incompatible with coarser boxes({refined_coarser_boxes}) and nest_buffer({nesting_buffer})")

            if box.ndim != ndim:
                print("box.ndim != ndim", box.ndim, ndim)
                raise ValueError(f"Box({box}) has incorrect dimensions for simulation")
            for l in boxm.refine(box, refinement_ratio).shape:
                if l < smallest_patch_size:
                    raise ValueError("Invalid box incompatible with smallest_patch_size")

    return refinement_boxes, len(refinement_boxes.items()) + 1


# ------------------------------------------------------------------------------

def check_patch_size(**kwargs):
    def get_max_ghosts():
        from ..core.gridlayout import GridLayout
        grid = GridLayout()
        return max(grid.nbrGhosts(kwargs["interp_order"], x) for x in ['primal','dual'])

    max_ghosts = get_max_ghosts()
    largest_patch_size = kwargs.get("largest_patch_size", None)
    smallest_patch_size = kwargs.get("smallest_patch_size", max_ghosts)

    cells = kwargs["cells"]

    if largest_patch_size is not None:

        if not all(size >= largest_patch_size for size in cells):
            raise ValueError("Error - largest_patch_size should be less than nbr of cells in all directions")

        if largest_patch_size <= 0:
            raise ValueError("Error - largest_patch_size cannot be <= 0")

    if not all(size >= smallest_patch_size for size in cells):
        raise ValueError("Error - smallest_patch_size should be less than nbr of cells in all directions")

    small_invalid_patch_size = (max_ghosts - 1)
    if smallest_patch_size <= small_invalid_patch_size:
        raise ValueError("Error - smallest_patch_size cannot be <= " + str(small_invalid_patch_size))

    if largest_patch_size is not None:
        if largest_patch_size < smallest_patch_size:
            raise ValueError("Error - largest_patch_size and smallest_patch_size are incompatible")

    return largest_patch_size, smallest_patch_size


# ------------------------------------------------------------------------------

# diag_options = {"format":"phareh5", "options": {"dir": "phare_ouputs/"}}
def check_diag_options(**kwargs):
    diag_options = kwargs.get("diag_options", None)
    formats = ["phareh5"]
    if diag_options is not None and "format" in diag_options:
        if diag_options["format"] not in formats:
            raise ValueError("Error - diag_options format is invalid")
        if "options" in diag_options and "dir" in diag_options["options"]:
            diag_dir = diag_options["options"]["dir"]
            if os.path.exists(diag_dir) and os.path.isfile(diag_dir):
                raise ValueError ("Error: Simulation diag_options dir exists as a file.")
            try:
                os.makedirs(diag_dir, exist_ok=True)
                if not os.path.exists(diag_dir):
                    raise ValueError ("1. Creation of the directory %s failed" % diag_dir)
            except FileExistsError:
                raise ValueError ("Creation of the directory %s failed" % diag_dir)
    return diag_options





def check_refinement(**kwargs):
    return kwargs.get("refinement", "boxes")



def check_nesting_buffer(**kwargs):
    nesting_buffer = kwargs.get('nesting_buffer', 0)

    if nesting_buffer < 0:
        raise ValueError(f"Error: nesting_buffer({nesting_buffer}) cannot be negative")

    smallest_patch_size = kwargs["smallest_patch_size"]
    largest_patch_size = kwargs["largest_patch_size"]

    if nesting_buffer > smallest_patch_size / 2:
        raise ValueError(f"Error: nesting_buffer({nesting_buffer})"
                          + f"cannot be larger than half the smallest_patch_size ({smallest_patch_size})")

    cells = np.asarray(kwargs["cells"])

    if largest_patch_size != None and (nesting_buffer > cells - kwargs["largest_patch_size"]).any():
      raise ValueError(f"Error: nesting_buffer({nesting_buffer})"
                        + f"incompatible with number of domain cells ({cells}) and largest_patch_size({largest_patch_size})")

    elif (nesting_buffer > cells).any():
            raise ValueError(f"Error: nesting_buffer({nesting_buffer}) incompatible with number of domain cells ({cells})")

    return nesting_buffer



def check_optional_keywords(**kwargs):
    extra = []

    if check_refinement(**kwargs) != "boxes":
        extra += ['max_nbr_levels']

    return extra



def check_resistivity(**kwargs):
    resistivity = kwargs.get("resistivity", 0.0)
    if resistivity < 0.0:
            raise ValueError(f"Error: resistivity should not be negative")

    return resistivity



def check_hyper_resistivity(**kwargs):
    hyper_resistivity = kwargs.get("hyper_resistivity", 0.0001)
    if hyper_resistivity < 0.0:
            raise ValueError(f"Error: hyper_resistivity should not be negative")

    return hyper_resistivity



# ------------------------------------------------------------------------------

def checker(func):
    def wrapper(simulation_object, **kwargs):
        accepted_keywords = ['domain_size', 'cells', 'dl', 'particle_pusher', 'final_time',
                             'time_step', 'time_step_nbr', 'layout', 'interp_order', 'origin',
                             'boundary_types', 'refined_particle_nbr', 'path', 'nesting_buffer',
                             'diag_export_format', 'refinement_boxes', 'refinement', 'init_time',
                             'smallest_patch_size', 'largest_patch_size', "diag_options",
                             'resistivity', 'hyper_resistivity' ]

        accepted_keywords += check_optional_keywords(**kwargs)

        wrong_kwds = phare_utilities.not_in_keywords_list(accepted_keywords, **kwargs)
        if len(wrong_kwds) > 0:
            raise ValueError("Error: invalid arguments - " + " ".join(wrong_kwds))

        dl, cells = check_domain(**kwargs)

        kwargs["refinement_ratio"] = 2

        kwargs["dl"] = dl
        kwargs["cells"] =  cells
        kwargs["refinement_ratio"] = 2

        kwargs["init_time"] = kwargs.get('init_time', 0)

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

        kwargs["boundary_types"] = check_boundaries(ndim, **kwargs)
        kwargs["origin"] = check_origin(ndim, **kwargs)

        kwargs["refined_particle_nbr"] = check_refined_particle_nbr(ndim, **kwargs)
        kwargs["diag_export_format"] = kwargs.get('diag_export_format', 'hdf5')
        assert kwargs["diag_export_format"] in ["hdf5"] # only hdf5 supported for now

        largest, smallest = check_patch_size(**kwargs)
        kwargs["smallest_patch_size"] = smallest
        kwargs["largest_patch_size"] = largest

        kwargs["nesting_buffer"] = check_nesting_buffer(**kwargs)

        kwargs["refinement"] = check_refinement(**kwargs)
        if kwargs["refinement"] == "boxes":
            kwargs["refinement_boxes"], kwargs["max_nbr_levels"] = check_refinement_boxes(ndim, **kwargs)
        else:
            kwargs["max_nbr_levels"] = kwargs.get('max_nbr_levels', None)
            assert kwargs["max_nbr_levels"] != None # this needs setting otherwise
            kwargs["refinement_boxes"] = None

        kwargs["resistivity"] = check_resistivity(**kwargs)

        kwargs["hyper_resistivity"] = check_hyper_resistivity(**kwargs)

        return func(simulation_object, **kwargs)

    return wrapper


# ------------------------------------------------------------------------------


class Simulation(object):
    """
    1D run example: Simulation(time_step_nbr = 100, boundary_types="periodic", cells=80)
    2D run example: Simulation(time_step_nbr = 100, boundary_types=("periodic","periodic"), cells=(80,20))
    3D run example: Simulation(time_step_nbr = 100, boundary_types=("periodic","periodic","periodic"), cells=(80,20,40))

    optional parameters:
    -------------------

    dl                   : grid spacing dx, (dx,dy) or (dx,dy,dz) in 1D, 2D or 3D
    domain_size          : size of the physical domain Lx, (Lx,Ly), (Lx,Ly,Lz) in 1D, 2D or 3D
    cells                : number of cells nx or (nx, ny) or (nx, ny, nz) in 1, 2 and 3D.
    final_time           : final simulation time. Must be set if 'time_step' is not
    time_step            : simulation time step. Must be specified if 'final_time' is not
    interp_order         : order of the particle/mesh interpolation. Either 1, 2 or 3 (default=1)
    layout               : layout of the physical quantities on the mesh (default = "yee")
    origin               : origin of the physical domain, (default (0,0,0) in 3D)
    refined_particle_nbr : number of refined particles for particle splitting ( TODO default hard-coded to 2)
    particle_pusher      : algo to push particles (default = "modifiedBoris")
    path                 : path for outputs (default : './')
    boundary_types       : type of boundary conditions (default is "periodic" for each direction)
    diag_export_format   : format of the output diagnostics (default= "phareh5")
    nesting_buffer       : [default=0] minimum gap in coarse cells from border between coarse and refined patch
    refinement_boxes     : [default=None] {"L0":{"B0":[(lox,loy,loz),(upx,upy,upz)],...,"Bi":[(),()]},..."Li":{B0:[(),()]}}
    smallest_patch_size  :
    largest_patch_size   :
    max_nbr_levels       : [default=1] max number of levels in the hierarchy if refinement_boxes != "boxes"
    init_time            : unused for now, will be time for restarts someday

    """

    @checker
    def __init__(self, **kwargs):

        if global_vars.sim is not None:
            raise RuntimeError("simulation is already created")

        global_vars.sim = self

        for k, v in kwargs.items():
            object.__setattr__(self, k, v)

        self.ndim = compute_dimension(self.cells)

        self.diagnostics = []
        self.model = None
        self.electrons = None

        # hard coded in C++ MultiPhysicsIntegrator::getMaxFinerLevelDt
        self.nSubcycles = 4
        self.stepDiff = 1/self.nSubcycles

        levelNumbers = list(range(self.max_nbr_levels))
        self.level_time_steps = [
          self.time_step * (self.stepDiff ** (ilvl)) for ilvl in levelNumbers
        ]
        self.level_step_nbr = [
          self.nSubcycles ** levelNumbers[ilvl] * self.time_step_nbr for ilvl in levelNumbers
        ]

    def final_time(self):
        return self.time_step * self.time_step_nbr

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





# ------------------------------------------------------------------------------

    def add_diagnostics(self, diag):
        if diag.name in [diag.name for diag in self.diagnostics]:
            raise ValueError("Error: diagnostics {} already registered".format(diag.name))

        # check whether the spatial extent of the diagnostics is valid, given the domain size
        if diag.extent() is not None and not self.within_simulation_extent(diag.extent()):
            raise RuntimeError("Error: invalid diagnostics spatial extent")

        self.diagnostics.append(diag)


# ------------------------------------------------------------------------------

    def count_diagnostics(self, type_name):
        return len([diag for diag in self.diagnostics if diag.type == type_name])


# ------------------------------------------------------------------------------

    def set_model(self, model):
        self.model = model


    def set_electrons(self, electrons):
        self.electrons = electrons
# ------------------------------------------------------------------------------
