
def basicSimulatorArgs(dim:int, interp:int, **kwargs):
    cells = [65 for i in range(dim)]
    if "cells" in kwargs:
        cells = kwargs["cells"]
    if not isinstance(cells, list):
        cells = [cells]
    dl = [1./v for v in cells]
    args = {
        "interp_order": interp,
        "smallest_patch_size":10,
        "largest_patch_size":64,
        "time_step_nbr":1000,
        "final_time":1.,
        "boundary_types":"periodic",
        "cells":cells,
        "dl":dl,
        "max_nbr_levels":2,
        "refinement_boxes" : {"L0":{"B0":[(10,),(50,)]}}
    }
    for k, v in kwargs.items():
        args[k] = v
    return args

def makeBasicModel():
    import phare.pharein as ph
    density = lambda x: 2.

    bx, by, bz = (lambda x: x if x != 0 else 1 for i in range(3))
    ex, ey, ez = (lambda x: x if x != 0 else 1 for i in range(3))
    vx, vy, vz = (lambda x: 1. for i in range(3))

    vthx, vthy, vthz = (lambda x: 1. for i in range(3))

    vvv = {
        "vbulkx":vx, "vbulky":vy, "vbulkz":vz,
        "vthx":vthx, "vthy":vthy, "vthz":vthz
    }

    ph.MaxwellianFluidModel(
        bx=bx, by=bx, bz=bx,
        ex=bx, ey=bx, ez=bx,
        protons={"charge":-1, "density":density, **vvv},
        alpha={"charge":-1, "density":density, **vvv}
    )
