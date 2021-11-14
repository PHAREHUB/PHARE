# default "job.py" for benchmark test suite

def _density(*xyz): return 1.0
def _bx(*xyz): return 1.0
def _by(*xyz): return 0.0
def _bz(*xyz): return 0.0
def _vx(*xyz): return 0.0
def _vy(*xyz): return 0.0
def _vz(*xyz): return 0.0
def _vthx(*xyz): return 0.30
def _vthy(*xyz): return 0.30
def _vthz(*xyz): return 0.30

def setup(**kwargs):
    from pyphare.core.phare_utilities import np_array_ify, NO_GUI
    import pyphare.pharein as ph
    NO_GUI()

    def getFn(key):
        return kwargs.get(key, globals()["_" + key])

    ndim = kwargs["ndim"]
    kwargs["cells"] = np_array_ify(kwargs["cells"], ndim)
    kwargs["dl"] = np_array_ify(kwargs["dl"], ndim)
    dl, cells = ph.simulation.check_domain(**kwargs)
    largest_patch_size, smallest_patch_size = ph.simulation.check_patch_size(**kwargs)
    _, time_step, final_time = ph.simulation.check_time(**kwargs)

    sim = ph.Simulation(
        interp_order=kwargs.get("interp_order", 1),
        smallest_patch_size=smallest_patch_size,
        largest_patch_size=largest_patch_size,
        time_step=time_step,
        final_time=final_time,
        boundary_types=[kwargs.get("boundary_types", "periodic")] * ndim,
        cells=cells,
        dl=np_array_ify(dl, ndim),
        diag_options={
            "format": "phareh5",
            "options": {"dir": kwargs.get("diag_dir", "."), "mode": "overwrite"},
        },
        **kwargs.get("kwargs", {})
    )
    # ph.simulation._print(sim)

    ph.MaxwellianFluidModel(
        bx=getFn("bx"),
        by=getFn("by"),
        bz=getFn("bz"),
        protons={
            "charge": kwargs.get("charge", 1),
            "density": getFn("density"),
            "nbr_part_per_cell": kwargs.get("ppc", 10),
            **{
                "vbulkx": getFn("vx"),
                "vbulky": getFn("vy"),
                "vbulkz": getFn("vz"),
                "vthx": getFn("vthx"),
                "vthy": getFn("vthy"),
                "vthz": getFn("vthz"),
            },
        },
    )
    ph.ElectronModel(closure="isothermal", Te=kwargs.get("Te", 0.12))
