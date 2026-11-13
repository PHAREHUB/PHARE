"""Small phareh5 files for tests, with the structure of real PHARE output.

Same groups, attributes and dtypes as src/diagnostic/detail/h5writer.hpp writes
(checked against build/sample_data/*/h5/*.h5): root attrs cell_width,
dimension, domain_box, interpOrder (hybrid only), layoutType, origin,
periodicities (+ pop_mass/pop_charge for population files), a py_attrs group,
``/t/<"%.10f">/pl<level>/p<rank>#<localid>/<dataset>`` with patch attrs lower,
upper, nbrCells, origin, mpi_rank, datasets C order ``[x, y, z]`` with ghost
nodes (attr ``ghosts``) at their true Yee positions.

Values come from a function of the physical coordinates of each stored value
(Yee positions, ghosts included). The default `linear_value` uses dyadic
coefficients so that, with dyadic mesh sizes, every stored value and every
2-point average is exact in float32: node values read back must be exactly
`linear_value` at the nodes.
"""

import numpy as np

DATASETS = {
    "EM_B": ("EM_B_x", "EM_B_y", "EM_B_z"),
    "EM_E": ("EM_E_x", "EM_E_y", "EM_E_z"),
    "ions_bulkVelocity": ("bulkVelocity_x", "bulkVelocity_y", "bulkVelocity_z"),
    "ions_charge_flux": ("charge_flux_x", "charge_flux_y", "charge_flux_z"),
    "ions_pop_protons_flux": ("flux_x", "flux_y", "flux_z"),
    "ions_charge_density": ("charge_density",),
    "ions_mass_density": ("mass_density",),
    "ions_pop_protons_density": ("density",),
    "ions_momentum_tensor": tuple(
        f"momentum_tensor_{c}" for c in ("xx", "xy", "xz", "yy", "yz", "zz")
    ),
    "mhd_rho": ("rho",),
    "mhd_V": ("V_x", "V_y", "V_z"),
}

# dual (cell-centred) directions per dataset name, from pyphare/core/gridlayout.py
_DUAL = {
    "EM_B_x": (False, True, True),
    "EM_B_y": (True, False, True),
    "EM_B_z": (True, True, False),
    "EM_E_x": (True, False, False),
    "EM_E_y": (False, True, False),
    "EM_E_z": (False, False, True),
}
_MHD_CELL = ("rho", "V_x", "V_y", "V_z", "P", "rhoV_x", "rhoV_y", "rhoV_z", "Etot")


def dual_of(dataset, dim, model="hybrid"):
    """Tuple of `dim` bools: True where `dataset` is dual (independent of the reader's table)."""
    if model == "mhd" and dataset in _MHD_CELL:
        return (True,) * dim
    return _DUAL.get(dataset, (False, False, False))[:dim]


def linear_coeffs(dataset):
    """(c0, (cx, cy, cz)) dyadic coefficients of the default linear field of `dataset`."""
    names = [d for ds in DATASETS.values() for d in ds]
    k = names.index(dataset) if dataset in names else 0
    c0 = 1.0 + (k % 7)
    coef = tuple(((k + a) % 4) - 1.5 for a in range(3))
    return c0, coef


def linear_value(dataset, coords):
    """Default field: c0 + sum_a c_a * coords[a]; coords = tuple of broadcastable arrays."""
    c0, coef = linear_coeffs(dataset)
    out = c0
    for a, x in enumerate(coords):
        out = out + coef[a] * x
    return out


def yee_coords(lower, nn, g, dual, dl):
    """1-D physical coordinates of the stored values of one dataset, ghosts included.

    lower : lower cell index per axis; nn : node counts (upper - lower + 2)
    g     : ghost width; dual : bools per axis; dl : mesh size per axis
    Primal index k is node lower + k - g, dual index k is cell lower + k - g
    (centre at + 0.5).
    """
    out = []
    for a in range(len(nn)):
        n = nn[a] - 1 + 2 * g if dual[a] else nn[a] + 2 * g
        idx = np.arange(n, dtype=np.float64) + lower[a] - g + (0.5 if dual[a] else 0.0)
        out.append(idx * dl[a])
    return out


def node_coords(lower, upper, dl):
    """1-D node coordinates per axis of a patch (nodes lower .. upper + 1)."""
    return [np.arange(lo, up + 2, dtype=np.float64) * d for lo, up, d in zip(lower, upper, dl)]


def write_phareh5(
    path,
    *,
    dim,
    patches,
    quantity="EM_B",
    datasets=None,
    func=linear_value,
    dl0=None,
    domain_cells=None,
    ghosts=2,
    dtype=np.float32,
    model="hybrid",
    nan_outer_ghost=False,
    chunked=False,
    pop_attrs=None,
):
    """Write a phareh5 file.

    patches       {time (float): {level: [(rank, localid, lower, upper), ...]}};
                  patch groups are created in the given order; a time whose
                  dict has no level 0 imitates a per-level dump
    quantity      key of DATASETS (or give `datasets` explicitly)
    func          func(dataset_name, coords) -> values, coords = tuple of `dim`
                  arrays (ij meshgrid of physical Yee positions incl. ghosts)
    dl0           level-0 mesh size per axis (default 0.5); level l uses dl0/2**l
    domain_cells  level-0 cells per axis (default 32)
    ghosts        ghost width written in every dataset (attr "ghosts")
    model         "hybrid" (writes interpOrder=1) or "mhd" (no interpOrder)
    nan_outer_ghost  put NaN in the outermost ghost layer of every dataset
    chunked       create chunked datasets (no single byte offset)
    pop_attrs     dict written as root attributes (e.g. pop_mass, pop_charge)
    """
    import h5py

    dl0 = np.full(dim, 0.5) if dl0 is None else np.asarray(dl0, dtype=np.float64)
    cells = np.full(dim, 32) if domain_cells is None else np.asarray(domain_cells)
    names = DATASETS[quantity] if datasets is None else tuple(datasets)
    with h5py.File(path, "w") as f:
        f.attrs["cell_width"] = dl0
        f.attrs["dimension"] = np.uint64(dim)
        f.attrs["domain_box"] = (cells - 1).astype(np.int32)
        if model == "hybrid":
            f.attrs["interpOrder"] = np.uint64(1)
        f.attrs["layoutType"] = "yee"
        f.attrs["origin"] = np.zeros(dim)
        f.attrs["periodicities"] = np.ones(dim, dtype=np.int32)
        for k, v in (pop_attrs or {}).items():
            f.attrs[k] = v
        pa = f.create_group("py_attrs")
        for k in ("git_hash", "highfive_version", "pybind_version", "samrai_version"):
            pa.attrs[k] = "test"
        pa.attrs["serialized_simulation"] = ""
        tg = f.create_group("t")
        for time, levels in patches.items():
            ttg = tg.create_group(f"{time:.10f}")
            for lvl, plist in levels.items():
                lg = ttg.create_group(f"pl{lvl}")
                dl = dl0 / 2.0**lvl
                for rank, localid, lower, upper in plist:
                    lower = np.asarray(lower, dtype=np.int64)
                    upper = np.asarray(upper, dtype=np.int64)
                    nn = upper - lower + 2
                    pg = lg.create_group(f"p{rank}#{localid}")
                    pg.attrs["lower"] = lower.astype(np.int32)
                    pg.attrs["mpi_rank"] = np.uint64(rank)
                    pg.attrs["nbrCells"] = (upper - lower + 1).astype(np.uint32)
                    pg.attrs["origin"] = lower * dl
                    pg.attrs["upper"] = upper.astype(np.int32)
                    for dname in names:
                        dual = dual_of(dname, dim, model)
                        axes = yee_coords(lower, nn, ghosts, dual, dl)
                        coords = np.meshgrid(*axes, indexing="ij")
                        data = np.asarray(func(dname, tuple(coords)), dtype=np.float64)
                        data = np.broadcast_to(data, coords[0].shape).astype(dtype)
                        if nan_outer_ghost and ghosts > 0:
                            for a in range(dim):
                                sl = [slice(None)] * dim
                                for edge in (0, -1):
                                    sl[a] = edge
                                    data[tuple(sl)] = np.nan
                        chunks = tuple(max(1, s // 2) for s in data.shape)
                        kw = {"chunks": chunks} if chunked else {}
                        ds = pg.create_dataset(dname, data=data, **kw)
                        ds.attrs["ghosts"] = np.uint64(ghosts)
    return path


def tiling(dim, cells, size, rank_of=None, first_id=0):
    """Patches (rank, localid, lower, upper) tiling `cells` (per axis) by blocks of `size`.

    rank_of(i) gives the rank of the i-th block (default 0); ids increase from
    `first_id` in C order of the blocks.
    """
    cells = np.broadcast_to(np.asarray(cells), (dim,))
    size = np.broadcast_to(np.asarray(size), (dim,))
    starts = [np.arange(0, c, s) for c, s in zip(cells, size)]
    grid = np.stack(np.meshgrid(*starts, indexing="ij"), axis=-1).reshape(-1, dim)
    out = []
    for i, lo in enumerate(grid):
        up = np.minimum(lo + size, cells) - 1
        rank = 0 if rank_of is None else rank_of(i)
        out.append((rank, first_id + i, tuple(lo.tolist()), tuple(up.tolist())))
    return out
