# default "job.py" for benchmark test suite

def _density(*xyz):
    return 1.0


def _bx(*xyz):
    return 1.0


def _by(*xyz):
    return 0.0


def _bz(*xyz):
    return 0.0


def _vx(*xyz):
    return 0.0


def _vy(*xyz):
    return 0.0


def _vz(*xyz):
    return 0.0


def _vthx(*xyz):
    return 0.30


def _vthy(*xyz):
    return 0.30


def _vthz(*xyz):
    return 0.30


def setup(**kwargs):
    from pyphare.pharein import Simulation, MaxwellianFluidModel
    from pyphare.pharein import ElectronModel, ElectromagDiagnostics, FluidDiagnostics
    import matplotlib as mpl

    mpl.use("Agg")

    def getFn(key):
        return kwargs.get(key, globals()["_" + key])

    density = getFn("density")
    bx, by, bz = getFn("bx"), getFn("by"), getFn("bz")
    vx, vy, vz = getFn("vx"), getFn("vy"), getFn("vz")
    vthx, vthy, vthz = getFn("vthx"), getFn("vthy"), getFn("vthz")


    ndim = kwargs.get("ndim", 1)
    interp = kwargs.get("ndim", 1)
    smallest_patch_size = kwargs.get("smallest_patch_size", 20)
    largest_patch_size = kwargs.get("largest_patch_size", 20)
    time_step_nbr = kwargs.get("time_step_nbr", 100)
    final_time = kwargs.get("final_time", 0.1)
    cells = [kwargs.get("cells", 100)] * ndim
    dl = [kwargs.get("dl", 0.2)] * ndim
    Te = kwargs.get("Te", 0.12)
    charge = kwargs.get("charge", 1)
    ppc = kwargs.get("ppc", 10)
    boundary_types = [kwargs.get("boundary_types", "periodic")] * ndim
    diag_dir = kwargs.get("diag_dir", ".")

    Simulation(
        interp_order=interp,
        smallest_patch_size=smallest_patch_size,
        largest_patch_size=largest_patch_size,
        time_step_nbr=time_step_nbr,
        final_time=final_time,
        boundary_types=boundary_types,
        cells=cells,
        dl=dl,
        diag_options={
            "format": "phareh5",
            "options": {"dir": diag_dir, "mode": "overwrite"},
        },
    )
    MaxwellianFluidModel(
        bx=bx,
        by=by,
        bz=bz,
        protons={
            "charge": charge,
            "density": density,
            "nbr_part_per_cell": ppc,
            **{
                "vbulkx": vx,
                "vbulky": vy,
                "vbulkz": vz,
                "vthx": vthx,
                "vthy": vthy,
                "vthz": vthz,
            },
        },
    )
    ElectronModel(closure="isothermal", Te=Te)
