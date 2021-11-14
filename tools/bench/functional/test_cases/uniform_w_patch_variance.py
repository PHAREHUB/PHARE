"""

  This test case is for assessing the impact of copy/stream for various patch size and parameter sets

  python3 -Ou tools/bench/functional/run_gen_plots.py --test_cases=uniform_w_patch_variance

"""
import os
from pathlib import Path
this_file_name, file_ext = os.path.splitext(os.path.basename(__file__))
gen_path = os.path.join(os.path.dirname(__file__), f'../generated/{this_file_name}')
#### DO NOT EDIT ABOVE ####

### test permutation section - minimized is best ###
time_step = .01
time_step_nbr = 1e2
dl = 0.25
cells = 1e2
ndims = [2]
mpirun_Ns = [1]
ppc_list = [10, 25, 50]
interps = [1, 2, 3]
patch_sizes = [25, 50, 100]
vth = { f"vth{xyz}" : lambda *xyz: .3 for xyz in "xyz"}

def generate(ndim, interp, ppc, mpirun_n, patch_size):
    """
      Params may include functions for the default population "protons"
         see: simulation_setup.py::setup for all available dict keys

         simulation_setup.setup doesn't even have to be used, any job.py style file is allowed
         A "params" dict must exist for exporting test case information
    """
    file_name = f"{int(ndim)}_{int(interp)}_{int(ppc)}_{int(mpirun_n)}_{int(patch_size)}"
    with open(os.path.join(gen_path, file_name + ".py"), "w") as out:
        out.write("""
import tools.bench.functional.test_cases.uniform_w_patch_variance as inputs # scary
params = {""" + f"""
    "mpirun_n"            : {mpirun_n},
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
""" + """
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
    for permutation in itertools.product(ndims, interps, ppc_list, mpirun_Ns, patch_sizes):
        generate(*permutation)

if __name__ == "__main__":
    generate_all()
