import os
from pathlib import Path

from pyphare.pharein.simulation import supported_dimensions

FILE_DIR = Path(__file__).resolve().parent
this_file_name, file_ext = os.path.splitext(os.path.basename(__file__))
gen_path = FILE_DIR / ".." / "generated" / this_file_name
#### DO NOT EDIT ABOVE ####


### test permutation section - minimized is best ###
dl = 0.2
time_step = 0.001
smallest_patch_size = 50
largest_patch_size = 50
time_step_nbr = 500

permutables = [
    ("ndim", supported_dimensions()),
    ("interp", [1, 2, 3]),
    ("cells", [50, 150, 300]),
    ("ppc", [50, 100, 200]),
]


def permutation_filename(ndim, interp, cells, ppc):
    return f"{ndim}_{interp}_{cells}_{ppc}.py"


def permutation_filepath(ndim, interp, cells, ppc):
    return str(gen_path / permutation_filename(ndim, interp, cells, ppc))


def generate(ndim, interp, cells, ppc):
    """
    Params may include functions for the default population "protons"
       see: simulation_setup.py::setup for all available dict keys

       simulation_setup.setup doesn't even have to be used, any job.py style file is allowed
       A "params" dict must exist for exporting test case information
    """
    filepath = permutation_filepath(ndim, interp, cells, ppc)
    with open(filepath, "w") as out:
        out.write(
            """
params = {"""
            + f"""
    "ndim"                : {ndim},
    "interp_order"        : {interp},
    "cells"               : {cells},
    "ppc"                 : {ppc},
    "time_step_nbr"       : {time_step_nbr},
    "dl"                  : {dl},
    "time_step"           : {time_step},
    "smallest_patch_size" : {smallest_patch_size},
    "largest_patch_size"  : {largest_patch_size}, """
            + """
}

import pyphare.pharein as ph
if ph.PHARE_EXE: # needed to allow params export without calling "job.py"
    from tools.bench.functional.simulation_setup import setup
    setup(**params) # basically a "job.py"

"""
        )


### following function is called during test_case generation ###
def generate_all(clean=True):
    gen_dir = Path(gen_path)
    if clean and os.path.exists(gen_path):
        import shutil

        shutil.rmtree(str(gen_dir))
    gen_dir.mkdir(parents=True, exist_ok=True)
    import itertools

    permutations = itertools.product(*[e[1] for e in permutables])
    for permutation in permutations:
        generate(*permutation)


if __name__ == "__main__":
    generate_all()
