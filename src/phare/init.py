#! /usr/bin/env python3

import phare.pyphare as pp
import sys, job

from phare.pharein.globals import sim as simulation

add = pp.add
addScalarFunction = getattr(pp, 'addScalarFunction{:d}'.format(simulation.dims)+'D')

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
add("simulation/time_step", float(simulation.time_step))
add("simulation/time_step_nbr",int(simulation.time_step_nbr))


add("simulation/AMR/max_nbr_levels", int(simulation.max_nbr_levels))
refinement_boxes = simulation.refinement_boxes



def as_paths(rb):
    add("simulation/AMR/refinement_boxes/nbr_levels/", int(len(rb.keys())))
    for level,boxes in rb.items():
        level_path = "simulation/AMR/refinement_boxes/"+level+"/"
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
                add(level_path+box_lower_path_y,format(lower[1]))
                add(level_path+box_upper_path_y,format(upper[1]))
                if (len(lower)==3):
                    box_lower_path_z = box + "/lower/z/"
                    box_upper_path_z = box + "/upper/z/"
                    add(level_path+box_lower_path_z, int(lower[2]))
                    add(level_path+box_upper_path_z, int(upper[2]))



if refinement_boxes is not None:
    as_paths(refinement_boxes)


add("simulation/algo/pusher/name", simulation.particle_pusher)

init_model = simulation.model
modelDict  = init_model.model_dict

add("simulation/ions/name", "ions")
add("simulation/ions/nbrPopulations", init_model.nbr_populations())


partinit = "particle_initializer"
for pop_index, pop in enumerate(init_model.populations):
    pop_path = "simulation/ions/pop"
    partinit_path = pop_path+"{:d}/".format(pop_index)+partinit+"/"
    d = modelDict[pop]
    add(pop_path+"{:d}/name".format(pop_index), pop)
    add(pop_path+"{:d}/mass".format(pop_index), float(d["mass"]))
    add(partinit_path+"name", "maxwellian")
    addScalarFunction(partinit_path+"density", d["density"])
    addScalarFunction(partinit_path+"bulk_velocity_x", d["vx"])
    addScalarFunction(partinit_path+"bulk_velocity_y", d["vy"])
    addScalarFunction(partinit_path+"bulk_velocity_z", d["vz"])
    addScalarFunction(partinit_path+"thermal_velocity_x",d["vthx"])
    addScalarFunction(partinit_path+"thermal_velocity_y",d["vthy"])
    addScalarFunction(partinit_path+"thermal_velocity_z",d["vthz"])
    add(partinit_path+"nbr_part_per_cell", int(d["nbrParticlesPerCell"]))
    add(partinit_path+"charge", float(d["charge"]))
    add(partinit_path+"basis", "cartesian")

add("simulation/electromag/name", "EM")
add("simulation/electromag/electric/name", "E")
elecinit_path = "simulation/electromag/electric/initializer/"
addScalarFunction(elecinit_path+"x_component", modelDict["ex"])
addScalarFunction(elecinit_path+"y_component", modelDict["ey"])
addScalarFunction(elecinit_path+"z_component", modelDict["ez"])

add("simulation/electromag/magnetic/name", "B")
maginit_path = "simulation/electromag/magnetic/initializer/"
addScalarFunction(maginit_path+"x_component", modelDict["bx"])
addScalarFunction(maginit_path+"y_component", modelDict["by"])
addScalarFunction(maginit_path+"z_component", modelDict["bz"])


diag_path = "simulation/diagnostics/"
for diag in simulation.diagnostics:
    categ_path = diag_path + diag.category + '/'
    name_path = categ_path + diag.name
    add(name_path + "/" + 'subtype/' , diag.diag_type)
    pp.add_size_t(name_path + "/" + 'compute_every/' , diag.compute_every)
    pp.add_size_t(name_path + "/" + 'write_every/' , diag.write_every)
    pp.add_size_t(name_path + "/" + 'start_iteration/' , diag.start_iteration)
    pp.add_size_t(name_path + "/" + 'last_iteration/' , diag.last_iteration)
add(diag_path + "filePath", "lol.5") # needs finishing

