import os

from core import phare_utilities
#from .globals import objects
from . import globals

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
            raise ValueError("Error: 'cells' and 'dl' must have the same dimension")

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

    return time_step_nbr, time_step


# ------------------------------------------------------------------------------


def check_interp_order(**kwargs):
    interp_order = kwargs.get('interp_order', 1)

    if interp_order not in [1, 2, 3, 4]:
        raise ValueError("Error: invalid interpolation order. Should be in [1,2,3,4]")

    return interp_order

# ------------------------------------------------------------------------------


def check_pusher(**kwargs):
    pusher = kwargs.get('particle_pusher', 'modified_boris')
    if pusher not in ['modified_boris']:
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


def check_boundaries(dims, **kwargs):
    valid_boundary_types = ("periodic",)
    boundary_types = kwargs.get('boundary_types', ['periodic'] * dims)
    phare_utilities.check_iterables(boundary_types)

    if phare_utilities.none_iterable(boundary_types):
        bc_length = 1
        if boundary_types not in valid_boundary_types:
            raise ValueError("Error: '{}' is not a valid boundary type".format(boundary_types))
        else:
            boundary_types = phare_utilities.listify(boundary_types)
    else:
        bc_length = len(boundary_types)
        for bc in boundary_types:
            if bc not in valid_boundary_types:
                raise ValueError("Error: '{}' is not a valid boundary type".format(bc))

    if bc_length != dims:
        raise ValueError("Error- boundary_types should have length {} and is of length {}".format(dims, bc_length))

    return boundary_types


# ------------------------------------------------------------------------------


def check_origin(dims, **kwargs):
    origin = kwargs.get("origin", [0.] * dims)
    return origin


# ------------------------------------------------------------------------------


def check_refinement_boxes(**kwargs):
    refinement_boxes = kwargs.get("refinement_boxes", None)
    if refinement_boxes is None:
        return refinement_boxes
    smallest_patch_size = kwargs.get("smallest_patch_size")
    for level,boxes in refinement_boxes.items():
        if (kwargs["max_nbr_levels"]) <= int(level[1:]) + 1:  # + 1 L0 index from 0 vs nbr_levels
            raise ValueError("Error - refinement_boxes, refined level larger than 'max_nbr_levels'")
        for box, points in boxes.items():
            for i in range(len(points[0])):
                if points[1][i] - points[0][i] + 1 < smallest_patch_size: # +1 as index from 0
                    raise ValueError("Invalid refineboxes, incompatible with smallest_patch_size")

    return refinement_boxes

# ------------------------------------------------------------------------------

def check_patch_size(**kwargs):
    def get_max_ghosts():
        from core.gridlayout import GridLayout
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
                if not os.path.exists(diag_dir):
                    os.makedirs(diag_dir, exist_ok=True)
                if not os.path.exists(diag_dir):
                    raise ValueError ("1. Creation of the directory %s failed" % diag_dir)
            except OSError:
                raise ValueError ("Creation of the directory %s failed" % diag_dir)
    return diag_options



# ------------------------------------------------------------------------------

def checker(func):
    def wrapper(simulation_object, **kwargs):
        accepted_keywords = ['domain_size', 'cells', 'dl', 'particle_pusher', 'final_time',
                             'time_step', 'time_step_nbr', 'layout', 'interp_order', 'origin',
                             'boundary_types', 'refined_particle_nbr', 'path',
                             'diag_export_format', 'max_nbr_levels', 'refinement_boxes',
                             'smallest_patch_size', 'largest_patch_size', "diag_options" ]

        wrong_kwds = phare_utilities.not_in_keywords_list(accepted_keywords, **kwargs)
        if len(wrong_kwds) > 0:
            raise ValueError("Error: invalid arguments - " + " ".join(wrong_kwds))

        dl, cells = check_domain(**kwargs)
        kwargs["dl"] = dl
        kwargs["cells"] =  cells

        time_step_nbr, time_step = check_time(**kwargs)
        kwargs["time_step_nbr"] = time_step_nbr
        kwargs["time_step"] = time_step

        kwargs["interp_order"] = check_interp_order(**kwargs)

        kwargs["particle_pusher"] = check_pusher(**kwargs)
        kwargs["layout"] = check_layout(**kwargs)
        kwargs["path"] = check_path(**kwargs)

        dims = compute_dimension(cells)
        kwargs["diag_options"] = check_diag_options(**kwargs)

        kwargs["boundary_types"] = check_boundaries(dims, **kwargs)
        kwargs["origin"] = check_origin(dims, **kwargs)

        kwargs["refined_particle_nbr"] = kwargs.get('refined_particle_nbr', 2)  # TODO change that default value
        kwargs["diag_export_format"] = kwargs.get('diag_export_format', 'ascii') #TODO add checker with valid formats

        largest, smallest = check_patch_size(**kwargs)
        kwargs["smallest_patch_size"] = smallest
        kwargs["largest_patch_size"] = largest

        kwargs["max_nbr_levels"] = kwargs.get('max_nbr_levels', 1)
        kwargs["refinement_boxes"] = check_refinement_boxes(**kwargs)

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
    interp_order         : order of the particle/mesh interpolation. Either 1, 2, 3 or 4 (default=1)
    layout               : layout of the physical quantities on the mesh (default = "yee")
    origin               : origin of the physical domain, (default (0,0,0) in 3D)
    refined_particle_nbr : number of refined particles for particle splitting ( TODO default hard-coded to 2)
    particle_pusher      : algo to push particles (default = "modifiedBoris")
    path                 : path for outputs (default : './')
    boundary_types       : type of boundary conditions (default is "periodic" for each direction)
    diag_export_format   : format of the output diagnostics (default= "phareh5")
    max_nbr_levels       : [default=1] max number of levels in the hierarchy
    refinement_boxes     : [default=None] {"L0":{"B0":[(lox,loy,loz),(upx,upy,upz)],...,"Bi":[(),()]},..."Li":{B0:[(),()]}}
    smallest_patch_size  :
    largest_patch_size   :

    """

    @checker
    def __init__(self, **kwargs):

        if globals.sim is not None:
            raise RuntimeError("simulation is already created")
        else:
            globals.sim = self
            # raise RuntimeError("simulation is already lol")

        for k, v in kwargs.items():
            object.__setattr__(self, k, v)

        self.dims = compute_dimension(self.cells)

        self.diagnostics = []
        self.model = None
        self.electrons = None

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
        else:
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
