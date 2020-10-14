

import numpy as np
import os, sys, subprocess

venv_path = os.environ.get('VIRTUAL_ENV')

if venv_path is not None:
    pythonexe = os.path.join(venv_path, "bin/python3")
    arg = "import sys; print(sys.path)"
    p = subprocess.run([pythonexe, "-c", arg], stdout=subprocess.PIPE)
    s = p.stdout.decode()
    s = s.replace('[','').replace(']','').replace('"',"").replace("\n","").replace("'",'')
    s = s.split(",")[1:]
    pythonpath = [ss.strip() for ss in s]
    sys.path = sys.path + pythonpath




from .uniform_model import UniformModel
from .maxwellian_fluid_model import MaxwellianFluidModel
from .electron_model import ElectronModel
from .diagnostics import FluidDiagnostics, ElectromagDiagnostics, ParticleDiagnostics
from .simulation import Simulation


def getSimulation():
    from .global_vars import sim
    return sim


def is_ndarray(x):
    return isinstance(x, np.ndarray)


def is_scalar(x):
    return not is_ndarray(x) and not isinstance(x, list)


# Wrap calls to user init functions to turn C++ vectors to ndarrays,
#  and returned ndarrays to C++ span
class fn_wrapper:

    def __init__(self, fn):
        self.fn = fn

    def __call__(self, *xyz):
        args = []

        for i, arg in enumerate(xyz):
            args.append(np.asarray(arg))

        ret = self.fn(*args)

        if isinstance(ret, list):
            ret = np.asarray(ret)

        if is_scalar(ret):
            ret = np.full(len(args[-1]), ret)

        from pybindlibs import cpp
        # convert numpy array to C++ SubSpan
        # couples vector init functions to C++
        return cpp.makePyArrayWrapper(ret)



def populateDict():

    from .global_vars import sim as simulation
    import pybindlibs.dictator as pp

    add = pp.add
    addInitFunction = getattr(pp, 'addInitFunction{:d}'.format(simulation.dims)+'D')

    add("simulation/name", "simulation_test")
    add("simulation/dimension", simulation.dims)
    add("simulation/boundary_types", simulation.boundary_types[0])

    if simulation.smallest_patch_size is not None:
        add("simulation/AMR/smallest_patch_size", simulation.smallest_patch_size)
    if simulation.largest_patch_size is not None:
        add("simulation/AMR/largest_patch_size", simulation.largest_patch_size)


    add("simulation/grid/layout_type", simulation.layout)
    add("simulation/grid/nbr_cells/x", int(simulation.cells[0]))
    add("simulation/grid/meshsize/x", float(simulation.dl[0]))
    add("simulation/grid/origin/x", float(simulation.origin[0]))

    if (simulation.dims>1):
        add("simulation/grid/nbr_cells/y", int(simulation.cells[1]))
        add("simulation/grid/meshsize/y", float(simulation.dl[1]))
        add("simulation/grid/origin/y", float(simulation.origin[1]))

        if (simulation.dims >2):
            add("simulation/grid/nbr_cells/z", int(simulation.cells[2]))
            add("simulation/grid/meshsize/z", float(simulation.dl[2]))
            add("simulation/grid/origin/z", float(simulation.origin[2]))



    add("simulation/interp_order", int(simulation.interp_order))
    add("simulation/refined_particle_nbr", int(simulation.refined_particle_nbr))
    add("simulation/time_step", float(simulation.time_step))
    add("simulation/time_step_nbr",int(simulation.time_step_nbr))


    add("simulation/AMR/max_nbr_levels", int(simulation.max_nbr_levels))
    refinement_boxes = simulation.refinement_boxes



    def as_paths(rb):
        add("simulation/AMR/refinement/boxes/nbr_levels/", int(len(rb.keys())))
        for level,boxes in rb.items():
            level_path = "simulation/AMR/refinement/boxes/"+level+"/"
            add(level_path + 'nbr_boxes/',int(len(boxes.keys())))
            for box,cells in boxes.items():
                lower = cells[0]
                upper = cells[1]
                box_lower_path_x = box + "/lower/x/"
                box_upper_path_x = box + "/upper/x/"
                add(level_path + box_lower_path_x, int(lower[0]))
                add(level_path + box_upper_path_x, int(upper[0]))
                if len(lower)>=2:
                    box_lower_path_y = box + "/lower/y/"
                    box_upper_path_y = box + "/upper/y/"
                    add(level_path+box_lower_path_y, int(lower[1]))
                    add(level_path+box_upper_path_y, int(upper[1]))
                    if (len(lower)==3):
                        box_lower_path_z = box + "/lower/z/"
                        box_upper_path_z = box + "/upper/z/"
                        add(level_path+box_lower_path_z, int(lower[2]))
                        add(level_path+box_upper_path_z, int(upper[2]))



    if refinement_boxes is not None and simulation.refinement =="boxes":
        as_paths(refinement_boxes)

    if simulation.refinement == "tagging":
        add("simulation/AMR/refinement/tagging/method","auto")


    add("simulation/algo/ion_updater/pusher/name", simulation.particle_pusher)

    init_model = simulation.model
    modelDict  = init_model.model_dict

    add("simulation/ions/nbrPopulations", init_model.nbr_populations())


    partinit = "particle_initializer"
    for pop_index, pop in enumerate(init_model.populations):
        pop_path = "simulation/ions/pop"
        partinit_path = pop_path+"{:d}/".format(pop_index)+partinit+"/"
        d = modelDict[pop]
        add(pop_path+"{:d}/name".format(pop_index), pop)
        add(pop_path+"{:d}/mass".format(pop_index), float(d["mass"]))
        add(partinit_path+"name", "maxwellian")

        addInitFunction(partinit_path+"density", fn_wrapper(d["density"]))
        addInitFunction(partinit_path+"bulk_velocity_x", fn_wrapper(d["vx"]))
        addInitFunction(partinit_path+"bulk_velocity_y", fn_wrapper(d["vy"]))
        addInitFunction(partinit_path+"bulk_velocity_z", fn_wrapper(d["vz"]))
        addInitFunction(partinit_path+"thermal_velocity_x",fn_wrapper(d["vthx"]))
        addInitFunction(partinit_path+"thermal_velocity_y",fn_wrapper(d["vthy"]))
        addInitFunction(partinit_path+"thermal_velocity_z",fn_wrapper(d["vthz"]))
        add(partinit_path+"nbr_part_per_cell", int(d["nbrParticlesPerCell"]))
        add(partinit_path+"charge", float(d["charge"]))
        add(partinit_path+"basis", "cartesian")
        if "init" in d and "seed" in d["init"]:
            pp.add_optional_size_t(partinit_path+"init/seed", d["init"]["seed"])

    add("simulation/electromag/name", "EM")
    add("simulation/electromag/electric/name", "E")


    add("simulation/electromag/magnetic/name", "B")
    maginit_path = "simulation/electromag/magnetic/initializer/"
    addInitFunction(maginit_path+"x_component", fn_wrapper(modelDict["bx"]))
    addInitFunction(maginit_path+"y_component", fn_wrapper(modelDict["by"]))
    addInitFunction(maginit_path+"z_component", fn_wrapper(modelDict["bz"]))

    diag_path = "simulation/diagnostics/"
    for diag in simulation.diagnostics:
        type_path = diag_path + diag.type + '/'
        name_path = type_path + diag.name
        add(name_path + "/" + 'type/' , diag.type)
        add(name_path + "/" + 'quantity/' , diag.quantity)
        pp.add_array_as_vector(name_path + "/" + "write_timestamps", diag.write_timestamps)
        pp.add_array_as_vector(name_path + "/" + "compute_timestamps", diag.compute_timestamps)


    if simulation.diag_options is not None and "options" in simulation.diag_options:
        add(diag_path + "filePath", simulation.diag_options["options"]["dir"])
    else:
        add(diag_path + "filePath", "phare_output")

    if (simulation.diag_options is not None and "options" in simulation.diag_options
        and "mode" in simulation.diag_options["options"]):
            add(diag_path + "mode", simulation.diag_options["options"]["mode"])


    #### adding electrons
    if simulation.electrons is None:
        raise RuntimeError("Error - no electrons registered to this Simulation")
    else:
        for item in simulation.electrons.dict_path():
            add("simulation/"+item[0], item[1])
