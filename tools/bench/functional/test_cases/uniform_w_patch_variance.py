"""
  This test case is for assessing the impact of copy/stream for various patch size and parameter sets
"""

import os
from pathlib import Path
from pyphare.pharein.simulation import supported_dimensions

FILE_DIR = Path(__file__).resolve().parent
this_file_name, file_ext = os.path.splitext(os.path.basename(__file__))
gen_path = FILE_DIR / ".." / "generated" / this_file_name
#### DO NOT EDIT ABOVE ####

### test permutation section - minimized is best ###
time_step = 0.001
time_step_nbr = 500
dl = 0.25
cells = 500

permutables = [
    ("ndim", supported_dimensions()),
    ("interp", [1, 2, 3]),
    ("ppc", [50, 100, 200]),
    ("patch_size", [25, 50, 100]),
]

vth = {f"vth{xyz}": lambda *xyz: 0.3 for xyz in "xyz"}


def permutation_filename(ndim, interp, ppc, patch_size):
    return f"{int(ndim)}_{int(interp)}_{int(ppc)}_{int(patch_size)}.py"


def permutation_filepath(ndim, interp, ppc, patch_size):
    return str(gen_path / permutation_filename(ndim, interp, ppc, patch_size))


def generate(ndim, interp, ppc, patch_size):
    """
    Params may include functions for the default population "protons"
       see: simulation_setup.py::setup for all available dict keys

       simulation_setup.setup doesn't even have to be used, any job.py style file is allowed
       A "params" dict must exist for exporting test case information
    """
    filepath = permutation_filepath(ndim, interp, ppc, patch_size)
    with open(filepath, "w") as out:
        out.write(
            """
import tools.bench.functional.test_cases.uniform_w_patch_variance as inputs # scary
params = {"""
            + f"""
    "ndim"                : {ndim},
    "interp_order"        : {interp},
    "ppc"                 : {ppc},
    "smallest_patch_size" : {patch_size},
    "largest_patch_size"  : {patch_size},
    "cells"               : inputs.cells,
    "time_step"           : inputs.time_step,
    "dl"                  : inputs.dl,
    "time_step_nbr"       : inputs.time_step_nbr,
    **inputs.vth,
"""
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
