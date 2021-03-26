import os
from pathlib import Path
from pyphare.pharein.simulation import supported_dimensions
this_file_name, file_ext = os.path.splitext(os.path.basename(__file__))
gen_path = os.path.join(os.path.dirname(__file__), f'../generated/{this_file_name}')
#### DO NOT EDIT ABOVE ####


### test permutation section - minimized is best ###
dl = 0.2
final_time = 0.1
smallest_patch_size = 20
largest_patch_size = 20
time_step_nbr = 1000
cells_list = [64, 128, 256]
ppc_list = [11, 111, 1111, 1111]
interps = [1, 2, 3]
ndims = supported_dimensions()


def generate(ndim, interp, cells, ppc):
    """
      Params may include functions for the default population "protons"
         see: simulation_setup.py::setup for all available dict keys

         simulation_setup.setup doesn't even have to be used, any job.py style file is allowed
         A "params" dict must exist for exporting test case information
    """
    file_name = f"{this_file_name}_{ndim}_{interp}_{cells}_{ppc}"
    with open(os.path.join(gen_path, file_name + ".py"), "w") as out:
        out.write("""
params = {""" + f"""
    "ndim"                : {ndim},
    "interp"              : {interp},
    "cells"               : {cells},
    "ppc"                 : {ppc},
    "time_step_nbr"       : {time_step_nbr},
    "dl"                  : {dl},
    "final_time"          : {final_time},
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
    for ndim in ndims:
        for interp in interps:
            for cells in cells_list:
                for ppc in ppc_list:
                    generate(ndim, interp, cells, ppc)


if __name__ == "__main__":
    generate_all()
