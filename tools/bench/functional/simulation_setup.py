# default "job.py" for benchmark test suite


def setup(**kwargs):
    from pyphare.core.phare_utilities import np_array_ify, NO_GUI
    import pyphare.pharein as ph

    NO_GUI()

    def density(*xyz):
        return 1.0

    def bx(*xyz):
        return 1.0

    def byz(*xyz):
        return 0.0

    def vxyz(*xyz):
        return 0.0

    def vthxyz(*xyz):
        return 0.30

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
    ph.MaxwellianFluidModel(
        bx=bx,
        by=byz,
        bz=byz,
        protons={
            "charge": kwargs.get("charge", 1),
            "density": density,
            "nbr_part_per_cell": kwargs.get("ppc", 10),
            **{
                "vbulkx": vxyz,
                "vbulky": vxyz,
                "vbulkz": vxyz,
                "vthx": vthxyz,
                "vthy": vthxyz,
                "vthz": vthxyz,
            },
        },
    )
    ph.ElectronModel(closure="isothermal", Te=kwargs.get("Te", 0.12))

    return sim
