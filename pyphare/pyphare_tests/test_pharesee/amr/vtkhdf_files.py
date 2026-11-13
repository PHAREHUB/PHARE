"""Write small .vtkhdf files with the exact structure of PHARE's vtkhdf writer.

Two generations of the writer are reproduced:

production  (files written before PR #1505, e.g. run0055Lx): float64 data,
            ``PointData/data`` chunks (1024, 1), ``AMRBox`` chunks (1024, 1),
            float64 ``Steps/Values``; old HDF5 file format (h5py default), so the
            chunk index is a version-1 B-tree; chunks allocated early, fill 0.
master      (current writer): float32 data, chunks (1 MiB // (ncomp * 4), ncomp)
            = (87381, 3) for vectors, ``AMRBox`` chunks (1024, 6), float32 times.

Layout reproduced from ``src/diagnostic/detail/vtkh5_type_writer.hpp`` and the
sample files of ``build/sample_data``: root attrs ``dimension`` (uint64),
``domain_box`` (int32, upper cell), ``interpOrder`` (uint64, hybrid only),
``periodicities`` (int32); ``/VTKHDF`` attrs Type, GridDescription (variable
length strings), Origin (float32 x3), Version (int32 [2, 2]); ``Level<l>`` attr
``Spacing`` (float32); each patch written ``copies`` times (4/2/1 in 1D/2D/3D),
nodes x fastest; per dump: Values appended, NSteps incremented, then per
level: AMRBox, NumberOfAMRBox, AMRBoxOffset, PointDataOffset, data.

Knobs produce the odd and broken files of DESIGN section 5.1 (see write_vtkhdf).
"""

import numpy as np

COPIES = {1: 4, 2: 2, 3: 1}


def default_fields(ncomp=3):
    """ncomp analytic functions f(t, level, *coords) -> values (float64).

    Smooth, distinct for every node, component, level and time, so that a
    misplaced value is detected.
    """

    def make(c):
        def f(t, level, *coords):
            x = coords[0]
            y = coords[1] if len(coords) > 1 else 0.0
            z = coords[2] if len(coords) > 2 else 0.0
            return (
                np.sin(1.7 * x + 0.3 + c) * np.cos(0.9 * y - 0.2 * c)
                + 0.37 * z
                + 0.01 * x
                + 0.001 * y
                + 10.0 * level
                + 100.0 * t
                + 1000.0 * c
            )

        return f

    return [make(c) for c in range(ncomp)]


def node_coords(lower, upper, dl):
    """Physical node coordinates of a patch, flattened x fastest: tuple of dim arrays."""
    lower = np.asarray(lower, dtype=np.int64)
    upper = np.asarray(upper, dtype=np.int64)
    axes = [np.arange(lo, up + 2) * d for lo, up, d in zip(lower, upper, dl)]
    grids = np.meshgrid(*axes, indexing="ij")
    return tuple(g.ravel(order="F") for g in grids)


def patch_values(fields, t, level, lower, upper, dl0, dtype):
    """(npts, ncomp) array of the analytic fields on the patch nodes, in `dtype`."""
    dl = np.asarray(dl0, dtype=np.float64) / 2**level
    coords = node_coords(lower, upper, dl)
    cols = [np.broadcast_to(np.asarray(f(t, level, *coords)), coords[0].shape) for f in fields]
    return np.stack(cols, axis=1).astype(dtype)


def _boxes(level_boxes, dim):
    """Normalise a level's boxes to (lower (n, dim), upper (n, dim)) int64."""
    if level_boxes is None or len(level_boxes) == 0:
        e = np.empty((0, dim), dtype=np.int64)
        return e, e
    arr = np.asarray(level_boxes, dtype=np.int64).reshape(len(level_boxes), 2, dim)
    return arr[:, 0, :], arr[:, 1, :]


def _create(group, name, shape, maxshape, dtype, chunks, compression=None):
    import h5py

    dcpl = h5py.h5p.create(h5py.h5p.DATASET_CREATE)
    if compression is None:
        dcpl.set_alloc_time(h5py.h5d.ALLOC_TIME_EARLY)
    return group.create_dataset(
        name,
        shape=shape,
        maxshape=maxshape,
        dtype=dtype,
        chunks=chunks,
        compression=compression,
        fillvalue=0,
        dcpl=dcpl,
    )


def _append(ds, values):
    values = np.asarray(values)
    n = ds.shape[0]
    ds.resize((n + len(values),) + ds.shape[1:])
    if len(values):
        ds[n : n + len(values)] = values


def write_vtkhdf(
    path,
    levels_boxes_per_step,
    fields=None,
    *,
    times=None,
    dl0=None,
    domain_cells=None,
    periodic=True,
    style="production",
    dtype=None,
    time_dtype=None,
    chunk_rows=None,
    copies=None,
    interp_order=1,
    attrs=None,
    libver=None,
    compression=None,
    offsets_dtype=np.int64,
    amrbox_cols=6,
    boundary_conditions=False,
    no_data=False,
    incomplete_last=None,
    emergency_dump=False,
    per_level_misalign=False,
    saturate_int32=False,
):
    """Write a .vtkhdf file and return its path.

    levels_boxes_per_step  list over steps of lists over levels of boxes; a
                    level's boxes are a sequence of (lower, upper) inclusive cell
                    boxes in that level's index space ([] = no patch at that step:
                    NumberOfAMRBox = 0). Levels missing at the end of a step are
                    padded with [].
    fields          list of ncomp callables f(t, level, *coords) (see
                    default_fields); default: 3 components
    times           one time per step (default 0, 0.5, 1, ...)
    dl0             level-0 mesh size per axis (default 0.4 everywhere)
    domain_cells    level-0 cells per axis (default: covers every level-0 box)
    style           "production" or "master" (module docstring); sets dtype,
                    time_dtype and chunk shapes unless given explicitly
    interp_order    int, or None for an MHD file (no interpOrder attribute)
    libver          passed to h5py.File (e.g. "latest": v2 headers, no v1 B-tree)
    compression     e.g. "gzip": compressed data chunks
    offsets_dtype   np.int64, or np.int32 as before PR #1287
    amrbox_cols     6 (unused columns 0) or 2*dim
    boundary_conditions  True: old attribute boundary_conditions=["periodic",...]
                    instead of periodicities
    no_data         True: only Steps/Values (what momentum_tensor files contain)
    incomplete_last None, "values" (one more time in Values, killed before any
                    level was written) or "level0" (one more step written for
                    level 0 only)
    emergency_dump  True: the last step's time is appended once per level, as
                    per-level emergency dumps do
    per_level_misalign  True: two more times written for levels >= 1 only
                    (fine_dump_lvl_max style): unsupported, FormatError expected
    saturate_int32  True: int32 offsets, the last PointDataOffset of level 0 is
                    2**31-1 (HDF5 clamps on overflow)
    """
    import h5py

    if style not in ("production", "master"):
        raise ValueError(style)
    production = style == "production"
    if fields is None:
        fields = default_fields(3)
    ncomp = len(fields)
    dtype = np.dtype(dtype or (np.float64 if production else np.float32))
    time_dtype = np.dtype(time_dtype or (np.float64 if production else np.float32))
    if saturate_int32:
        offsets_dtype = np.int32
    nsteps = len(levels_boxes_per_step)
    nlevels = max(len(s) for s in levels_boxes_per_step) if nsteps else 1
    first = next(b for s in levels_boxes_per_step for b in s if b is not None and len(b))
    dim = np.asarray(first[0][0]).reshape(-1).size
    copies = COPIES[dim] if copies is None else copies
    dl0 = np.full(dim, 0.4) if dl0 is None else np.asarray(dl0, dtype=np.float64)
    times = [0.5 * s for s in range(nsteps)] if times is None else list(times)
    if domain_cells is None:
        uppers = [_boxes(s[0], dim)[1] for s in levels_boxes_per_step if len(s[0])]
        domain_cells = np.max(np.concatenate(uppers), axis=0) + 1
    domain_cells = np.asarray(domain_cells, dtype=np.int64)
    if chunk_rows is None:
        chunk_rows = 1024 if production else max(1, (1 << 20) // (ncomp * dtype.itemsize))
    box_chunks = (1024, 1) if production else (1024, amrbox_cols)

    with h5py.File(path, "w", libver=libver) as f:
        f.attrs["dimension"] = np.uint64(dim)
        f.attrs["domain_box"] = (domain_cells - 1).astype(np.int32)
        if interp_order is not None:
            f.attrs["interpOrder"] = np.uint64(interp_order)
        per = np.broadcast_to(np.asarray(periodic, dtype=bool), (dim,))
        if boundary_conditions:
            f.attrs["boundary_conditions"] = [
                "periodic" if p else "physical" for p in per.tolist()
            ]
        else:
            f.attrs["periodicities"] = per.astype(np.int32)
        for k, v in (attrs or {}).items():
            f.attrs[k] = v
        root = f.create_group("VTKHDF")
        root.attrs["GridDescription"] = "XYZ"
        root.attrs["Origin"] = np.zeros(3, dtype=np.float32)
        root.attrs["Type"] = "OverlappingAMR"
        root.attrs["Version"] = np.array([2, 2], dtype=np.int32)
        steps = root.create_group("Steps")
        steps.attrs["NSteps"] = np.int32(0)
        values = _create(steps, "Values", (0,), (None,), time_dtype, (1024,))

        def add_time(t):
            _append(values, np.array([t], dtype=time_dtype))
            steps.attrs["NSteps"] = np.int32(int(steps.attrs["NSteps"]) + 1)

        if no_data:
            for t in times:
                add_time(t)
            return path

        groups = []
        for lvl in range(nlevels):
            g = root.create_group(f"Level{lvl}")
            g.attrs["Spacing"] = (dl0 / 2**lvl).astype(np.float32)
            g.create_group("CellData")
            g.create_group("FieldData")
            pd = g.create_group("PointData")
            data = _create(pd, "data", (0, ncomp), (None, ncomp), dtype, (chunk_rows, ncomp)
                           if not production else (chunk_rows, 1), compression)
            boxes = _create(g, "AMRBox", (0, amrbox_cols), (None, amrbox_cols), np.int32,
                            box_chunks)
            sl = steps.create_group(f"Level{lvl}")
            sl.create_group("CellDataOffset")
            sl.create_group("FieldDataOffset")
            tabs = {
                key: _create(sl, key, (0,), (None,), offsets_dtype, (1024,))
                for key in ("AMRBoxOffset", "NumberOfAMRBox", "PointDataOffset/data")
            }
            groups.append((data, boxes, tabs))

        def write_level(t, lvl, level_boxes):
            data, boxes, tabs = groups[lvl]
            lower, upper = _boxes(level_boxes, dim)
            _append(tabs["NumberOfAMRBox"], np.array([len(lower)], dtype=offsets_dtype))
            nbox = boxes.shape[0]
            _append(tabs["AMRBoxOffset"], np.array([nbox], dtype=offsets_dtype))
            flat = np.zeros((len(lower), amrbox_cols), dtype=np.int32)
            flat[:, 0 : 2 * dim : 2] = lower
            flat[:, 1 : 2 * dim : 2] = upper
            _append(boxes, flat)
            row0 = data.shape[0]
            _append(tabs["PointDataOffset/data"], np.array([row0], dtype=offsets_dtype))
            npts = np.prod(upper - lower + 2, axis=1) if len(lower) else np.zeros(0, int)
            nrows = int(npts.sum()) * copies
            # resize first (chunks allocated now, as the writer does), then write
            # in batches of ~64 MB to keep memory bounded on big files
            data.resize((row0 + nrows, ncomp))
            row, blocks, size = row0, [], 0
            for i, (lo, up) in enumerate(zip(lower, upper)):
                values = patch_values(fields, t, lvl, lo, up, dl0, dtype)
                blocks.append(np.tile(values, (copies, 1)))
                size += blocks[-1].nbytes
                if size > (64 << 20) or i == len(lower) - 1:
                    block = np.concatenate(blocks)
                    data[row : row + len(block)] = block
                    row += len(block)
                    blocks, size = [], 0

        def padded(step_boxes):
            return list(step_boxes) + [[]] * (nlevels - len(step_boxes))

        for s, (t, step_boxes) in enumerate(zip(times, levels_boxes_per_step)):
            last = s == nsteps - 1
            if emergency_dump and last:
                for lvl, lb in enumerate(padded(step_boxes)):
                    add_time(t)
                    write_level(t, lvl, lb)
                continue
            add_time(t)
            for lvl, lb in enumerate(padded(step_boxes)):
                write_level(t, lvl, lb)

        t_next = times[-1] + 1.0 if times else 1.0
        if incomplete_last == "values":
            add_time(t_next)
        elif incomplete_last == "level0":
            add_time(t_next)
            write_level(t_next, 0, padded(levels_boxes_per_step[-1])[0])
        elif incomplete_last is not None:
            raise ValueError(incomplete_last)
        if per_level_misalign:
            for k in range(2):
                tk = t_next + 0.25 * (k + 1)
                add_time(tk)
                for lvl in range(1, nlevels):
                    write_level(tk, lvl, padded(levels_boxes_per_step[-1])[lvl])
        if saturate_int32:
            ds = groups[0][2]["PointDataOffset/data"]
            ds[ds.shape[0] - 1] = np.iinfo(np.int32).max
    return path


def append_vtkhdf(path, levels_boxes_per_step, fields=None, *, times, dl0=None):
    """Append dumps to an existing file the way a running simulation does (h5py "a" mode:
    tables and data resized in place, new rows at the end). Same arguments as
    write_vtkhdf for the appended steps; dtypes and chunking are the file's."""
    import h5py

    fields = default_fields(3) if fields is None else fields
    with h5py.File(path, "a") as f:
        dim = int(f.attrs["dimension"])
        copies = COPIES[dim]
        dl0 = np.full(dim, 0.4) if dl0 is None else np.asarray(dl0, dtype=np.float64)
        root = f["VTKHDF"]
        steps = root["Steps"]
        nlevels = 0
        while f"Level{nlevels}" in root:
            nlevels += 1
        for t, step_boxes in zip(times, levels_boxes_per_step):
            values = steps["Values"]
            _append(values, np.array([t], dtype=values.dtype))
            steps.attrs["NSteps"] = np.int32(int(steps.attrs["NSteps"]) + 1)
            for lvl in range(nlevels):
                lb = step_boxes[lvl] if lvl < len(step_boxes) else []
                lower, upper = _boxes(lb, dim)
                data = root[f"Level{lvl}/PointData/data"]
                boxes = root[f"Level{lvl}/AMRBox"]
                st = steps[f"Level{lvl}"]
                odt = st["NumberOfAMRBox"].dtype
                _append(st["NumberOfAMRBox"], np.array([len(lower)], dtype=odt))
                _append(st["AMRBoxOffset"], np.array([boxes.shape[0]], dtype=odt))
                flat = np.zeros((len(lower), boxes.shape[1]), dtype=boxes.dtype)
                flat[:, 0 : 2 * dim : 2] = lower
                flat[:, 1 : 2 * dim : 2] = upper
                _append(boxes, flat)
                _append(st["PointDataOffset/data"], np.array([data.shape[0]], dtype=odt))
                rows = [np.tile(patch_values(fields, t, lvl, lo, up, dl0, data.dtype),
                                (copies, 1)) for lo, up in zip(lower, upper)]
                if rows:
                    _append(data, np.concatenate(rows))
    return path


def tiling(domain_cells, patch_cells, origin=None, stop=None):
    """Boxes tiling [origin, stop) with patches of `patch_cells` cells (edge patches smaller)."""
    domain_cells = np.asarray(domain_cells, dtype=np.int64)
    origin = np.zeros_like(domain_cells) if origin is None else np.asarray(origin, np.int64)
    stop = domain_cells if stop is None else np.asarray(stop, np.int64)
    patch_cells = np.broadcast_to(np.asarray(patch_cells, np.int64), domain_cells.shape)
    axes = [np.arange(o, s, p) for o, s, p in zip(origin, stop, patch_cells)]
    out = []
    for lo in np.stack(np.meshgrid(*axes, indexing="ij"), -1).reshape(-1, len(domain_cells)):
        up = np.minimum(lo + patch_cells, stop) - 1
        out.append((lo.tolist(), up.tolist()))
    return out


def naive_read(path, time_index, level):
    """Reference reader (plain h5py): {patch i: (npts, ncomp) array}, first copy only.

    time_index is the step index in the Steps tables (no dedup, no rounding).
    """
    import h5py

    with h5py.File(path, "r") as f:
        root = f["VTKHDF"]
        dim = int(f.attrs["dimension"])
        st = root[f"Steps/Level{level}"]
        nb = int(st["NumberOfAMRBox"][time_index])
        bo = int(st["AMRBoxOffset"][time_index])
        po = int(st["PointDataOffset/data"][time_index])
        boxes = root[f"Level{level}/AMRBox"][bo : bo + nb, : 2 * dim].astype(np.int64)
        data = root[f"Level{level}/PointData/data"]
        out = []
        row = po
        for b in boxes:
            shape = b[1::2] - b[0::2] + 2
            n = int(np.prod(shape))
            out.append((b[0::2], b[1::2], data[row : row + n, :]))
            row += n * COPIES[dim]
        return out
