import os
from pathlib import Path
from pyphare.pharein.simulation import supported_dimensions
this_file_name, file_ext = os.path.splitext(os.path.basename(__file__))
gen_path = os.path.join(os.path.dirname(__file__), f'../generated/{this_file_name}')
#### DO NOT EDIT ABOVE ####


### test permutation section - minimized is best ###
dl = 0.2
time_step = .01
time_step_nbr = 1e2
smallest_patch_size = 10
largest_patch_size = 50
cells_list = [10, 20, 50, 100, 200]
ppc_list = [10, 20, 50, 100, 200]
interps = [1, 2, 3]#
ndims = [2] #supported_dimensions()
threads = [10]


def generate(ndim, interp, cells, ppc, threads):
    """
      Params may include functions for the default population "protons"
         see: simulation_setup.py::setup for all available dict keys

         simulation_setup.setup doesn't even have to be used, any job.py style file is allowed
         A "params" dict must exist for exporting test case information
    """
    file_name = f"{this_file_name}_{ndim}_{interp}_{cells}_{int(ppc)}_{threads}"
    with open(os.path.join(gen_path, file_name + ".py"), "w") as out:
        out.write("""
params = {""" + f"""
    "ndim"                : {ndim},
    "interp_order"        : {interp},
    "cells"               : {cells},
    "ppc"                 : {ppc},
    "time_step"           : {time_step},
    "time_step_nbr"       : {time_step_nbr},
    "dl"                  : {dl},
    "threads"             : {threads},
    "smallest_patch_size" : {smallest_patch_size},
    "largest_patch_size"  : {largest_patch_size}, """ + """
}

import pyphare.pharein as ph
if ph.PHARE_EXE: # needed to allow params export without calling "job.py"
    from tools.bench.functional.simulation_setup import setup
    setup(**params) # basically a "job.py"

"""     )

### following function is called during test_case generation ###
def generate_all(clean=True):
    gen_dir = Path(gen_path)
    if clean and os.path.exists(gen_path):
        import shutil
        shutil.rmtree(str(gen_dir))
    gen_dir.mkdir(parents=True, exist_ok=True)
    import itertools
    for permutation in itertools.product(ndims, interps, cells_list, ppc_list, threads):
        generate(*permutation)


if __name__ == "__main__":
    generate_all()
