import os
from pathlib import Path

from pyphare.pharein.simulation import supported_dimensions

this_file_name, file_ext = os.path.splitext(os.path.basename(__file__))
gen_path = os.path.join(os.path.dirname(__file__), f'../generated/{this_file_name}')

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
if ph.PHARE_EXE:
    from tools.bench.functional.simulation_setup import setup
    setup(**params)

"""     )


def generate_all(clean=False):
    gen_dir = Path(gen_path)
    if clean:
        gen_dir.unlink()
    if not os.path.exists(gen_path):
        gen_dir.mkdir(parents=True, exist_ok=True)
        for ndim in ndims:
            for interp in interps:
                for cells in cells_list:
                    for ppc in ppc_list:
                        generate(ndim, interp, cells, ppc)


if __name__ == "__main__":
    generate_all()
