"""Point probes and line cuts of an AMR snapshot.

A point takes its value from the FINEST level that has a patch whose closed node box
contains it, by multilinear interpolation on the 2^d nodes of the cell of that patch that
contains it. A node shared by several patches of that level takes its value from ONE of
them: the finest level, then the lowest store slot (the weights of the nodes off the face
are exactly 0, so only the copies on the face matter). Copies of a shared node are equal in
single-rank output, but in multi-rank PHARE output copies held by patches of different MPI
ranks can differ at the 1e-4 relative level (measured 9e-5 .. 1.8e-4 on the 3-rank sample
run): the probe then returns one of them, not their mean. ``field.check_shared_nodes()``
measures that spread per level.

Values are returned in float64 (``np.result_type(field dtype, float64)``): interpolated
values are new numbers, they are not rounded back to the file's float32.
"""

from __future__ import annotations

import itertools

import numpy as np

from ..boxes import AXES, fortran_strides
from .finest import _flat, _names

# --------------------------------------------------------------------------- points


def _snap_int(q):
    """Snap level-index coordinates within 1e-9 (relative) of an integer onto it.

    Why: 0.6 / 0.2 == 2.9999999999999996 must be node 3, so that a point given at a node
    coordinate lands exactly on it (weights exactly 0 and 1).
    """
    r = np.round(q)
    return np.where(np.abs(q - r) <= 1e-9 * np.maximum(1.0, np.abs(q)), r, q)


def _prepare(hier, pts):
    """Validate points ``(N, dim)`` float64: wrap periodic axes, check domain and extent."""
    pts = np.array(pts, dtype=np.float64, copy=True)
    dim = hier.dim
    if pts.ndim != 2 or pts.shape[1] != dim:
        raise ValueError(f"points must be an array (N, {dim}) of coordinates, got {pts.shape}")
    bad = ~np.isfinite(pts).all(axis=1)
    if bad.any():
        raise ValueError(f"point {pts[np.argmax(bad)].tolist()} is not finite")
    domain = np.asarray(hier.domain, np.float64)
    for a in range(dim):
        out = (pts[:, a] < 0) | (pts[:, a] > domain[a])
        if not out.any():
            continue
        # DESIGN-NOTE: only points strictly outside [0, L] are wrapped; x == L stays (the same
        # node as x == 0 in periodic data).
        if hier.periodic[a]:
            pts[out, a] = np.mod(pts[out, a], domain[a])
        else:
            i = int(np.argmax(out))
            raise ValueError(
                f"point {pts[i].tolist()} is outside the domain [0, {domain[a]}] along "
                f"{AXES[a]}, which is not periodic; give points inside the domain"
            )
    if hier.extent is not None:
        lo, hi = hier.extent.clipped(domain)
        out = ((pts < lo) | (pts > hi)).any(axis=1)
        if out.any():
            i = int(np.argmax(out))
            raise ValueError(
                f"point {pts[i].tolist()} is outside the extent of this field "
                f"{hier.extent}; load a larger box or probe inside it"
            )
    return pts


def _cells_around(xi, ncells):
    """Cell range ``[ceil(xi) - 1, floor(xi)]`` clipped to the level: a patch's closed node
    box ``[lower, upper + 1]`` contains ``xi`` iff its cell box meets this range."""
    lo = np.clip(np.ceil(xi).astype(np.int64) - 1, 0, ncells - 1)
    hi = np.clip(np.floor(xi).astype(np.int64), 0, ncells - 1)
    return lo, hi


def locate(hier, pts):
    """Finest level and patch holding every point.

    pts : (N, dim) float64 physical coordinates, inside the domain (see ``_prepare``).
    Returns ``(ilev, slot, xi)``: ``ilev`` (N,) int64 index into ``hier.levels`` (-1 if no
    level holds the point), ``slot`` (N,) store slot in that level, ``xi`` (N, dim) float64
    the point in that level's node index units. Vectorised (BinIndex pair queries), finest
    level first, only the points not found yet go to the coarser level. No data is read.
    """
    n = len(pts)
    ilev = np.full(n, -1, np.int64)
    slot = np.full(n, -1, np.int64)
    xi_out = np.zeros_like(pts)
    todo = np.arange(n, dtype=np.int64)
    for k in range(len(hier.levels) - 1, -1, -1):
        if len(todo) == 0:
            break
        geom = hier.levels[k]
        if geom.npatch == 0:
            continue
        xi = _snap_int(pts[todo] / geom.dl)
        qlo, qhi = _cells_around(xi, geom.domain_cells)
        pairs = geom.bin_index().pairs(qlo, qhi)
        # pairs are sorted by query row: the first pair of each row is its patch
        rows, first = np.unique(pairs.a, return_index=True)
        found = todo[rows]
        ilev[found] = k
        slot[found] = pairs.b[first]
        xi_out[found] = xi[rows]
        todo = np.setdiff1d(todo, found, assume_unique=True)
    return ilev, slot, xi_out


def _interp(geom, flat, slots, xi, out_dtype):
    """Multilinear interpolation of ``flat`` at ``xi`` (level node units) inside ``slots``."""
    lower = geom.lower[slots]
    shape = geom.shape[slots]
    stride = fortran_strides(shape)
    u = xi - lower
    i0 = np.clip(np.floor(u).astype(np.int64), 0, shape - 2)
    t = u - i0
    base = geom.offset[slots] + (i0 * stride).sum(axis=1)
    val = np.zeros(len(slots), dtype=out_dtype)
    for corner in itertools.product((0, 1), repeat=geom.dim):
        c = np.array(corner, dtype=np.int64)
        w = np.prod(np.where(c == 1, t, 1.0 - t), axis=1)
        val += w * flat[base + (stride * c).sum(axis=1)]
    return val


def sample(snap, pts, names=None):
    """Values of the snapshot at points ``pts`` (N, dim), float64 dict ``name -> (N,)``.

    Reads (``ensure``) only the patches holding a point, one call per level.
    Raises ValueError for points outside the domain (non-periodic axis) or the extent.
    """
    hier = snap.hier
    sel = _names(snap, names)
    pts = _prepare(hier, pts)
    ilev, slot, xi = locate(hier, pts)
    if (ilev < 0).any():  # cannot happen when level 0 tiles the domain
        i = int(np.argmax(ilev < 0))
        raise ValueError(f"no patch holds point {pts[i].tolist()}")
    # DESIGN-NOTE: probe values are float64 (interpolated numbers, not data), see module doc.
    out_dtype = np.result_type(snap.levels[0].dtype, np.float64)
    out = {name: np.empty(len(pts), dtype=out_dtype) for name in sel}
    for k, (geom, ld) in enumerate(zip(hier.levels, snap.levels)):
        rows = np.nonzero(ilev == k)[0]
        if len(rows) == 0:
            continue
        ld.ensure(np.unique(slot[rows]), sel)
        for name in sel:
            out[name][rows] = _interp(geom, _flat(ld, name), slot[rows], xi[rows], out_dtype)
    return out


# --------------------------------------------------------------------------- lines


def line_coords(hier, fixed, free):
    """Coordinates along axis ``free`` of the finest nodes available on a line.

    fixed : dict ``axis -> coordinate`` for the other axes (already wrapped/validated).
    Every level contributes the nodes of its patches whose closed node box contains the
    line, except nodes strictly inside an interval covered by a finer level's patches. The
    result is sorted, unique (computed on integers of the finest level), restricted to the
    extent, float64.
    """
    top = hier.levels[-1].level
    keep_parts = []
    finer = []  # (A, B) node intervals on the finest grid, covered by finer levels
    for geom in reversed(hier.levels):
        r = 2 ** (top - geom.level)
        lo = np.zeros(geom.dim, np.int64)
        hi = geom.domain_cells - 1
        for a, v in fixed.items():
            xi = _snap_int(np.array([v / geom.dl[a]]))
            clo, chi = _cells_around(xi, geom.domain_cells[a])
            lo[a], hi[a] = clo[0], chi[0]
        slots = geom.patches_touching(lo, hi)
        if len(slots) == 0:
            continue
        a0 = geom.lower[slots, free] * r
        b0 = (geom.upper[slots, free] + 1) * r
        # every node of every crossing patch, on the finest integer grid
        counts = geom.shape[slots, free]
        starts = np.repeat(geom.lower[slots, free], counts)
        local = np.arange(int(counts.sum()), dtype=np.int64) - np.repeat(
            np.cumsum(counts) - counts, counts
        )
        nodes = (starts + local) * r
        if finer:
            A = np.concatenate([f[0] for f in finer])
            B = np.concatenate([f[1] for f in finer])
            order = np.argsort(A, kind="stable")
            A, Bmax = A[order], np.maximum.accumulate(B[order])
            j = np.searchsorted(A, nodes, side="left") - 1  # last interval with A < node
            inside = (j >= 0) & (Bmax[np.maximum(j, 0)] > nodes)
            nodes = nodes[~inside]
        keep_parts.append(nodes)
        finer.append((a0, b0))
    nodes = np.unique(np.concatenate(keep_parts)) if keep_parts else np.zeros(0, np.int64)
    dl_top = hier.levels[-1].dl[free]
    coords = nodes * dl_top
    if hier.extent is not None:
        lo, hi = hier.extent.clipped(hier.domain)
        coords = coords[(coords >= lo[free]) & (coords <= hi[free])]
    return coords


# --------------------------------------------------------------------------- public


def probe(snap, x=None, y=None, z=None, points=None):
    """Sample a snapshot at points or along a line.

    Points (multilinear interpolation in the finest level holding the point):
      ``probe(snap, points=P)``    P (N, dim) coordinates -> values (N,)
      ``probe(snap, x=xs, y=ys)``  one coordinate per axis, scalars or equal-length arrays ->
                                   values with the broadcast shape of the inputs
    Line (one free axis):
      2D ``probe(snap, y=800.)`` (or ``x=``) -> ``(coords, values)``: ``coords`` the node
      coordinates of the finest data along the line (each segment at its own level's
      spacing, sorted, no duplicates: NOT uniform), values interpolated across the fixed
      coordinate; 3D needs two fixed coordinates.
    Values are float64 arrays for a one-component snapshot, else a dict ``name -> array``.
    Coordinates are physical (same units as the node coordinates); a point outside the
    domain is wrapped on a periodic axis. Reads only the patches holding the points.
    Raises ValueError (point outside the domain on a non-periodic axis, or outside the
    field's extent; the message names the point), TypeError (bad combination of
    arguments, e.g. nothing fixed in 1D), NotImplementedError (a plane in 3D: use finest()).
    """
    hier = snap.hier
    dim = hier.dim
    single = snap.names[0] if len(snap.names) == 1 else None

    def pack(vals):
        return vals[single] if single is not None else vals

    given = {a: v for a, v in enumerate((x, y, z)) if v is not None}
    if any(a >= dim for a in given):
        extra = [AXES[a] for a in given if a >= dim]
        raise TypeError(f"{extra} given but the data is {dim}D (axes {AXES[:dim]})")
    if points is not None:
        if given:
            raise TypeError("give either points= or coordinates x=, y=, z=, not both")
        pts = np.asarray(points, dtype=np.float64)
        if dim == 1 and pts.ndim == 1:
            pts = pts[:, None]
        return pack(sample(snap, pts))
    # DESIGN-NOTE: x=, y=, z= scalars give a 0-d result; arrays broadcast together and the
    # result has their broadcast shape.
    if len(given) == dim:
        arrays = np.broadcast_arrays(*[np.asarray(given[a], np.float64) for a in range(dim)])
        shape = arrays[0].shape
        pts = np.stack([a.ravel() for a in arrays], axis=1)
        vals = sample(snap, pts)
        return pack({k: v.reshape(shape) for k, v in vals.items()})
    if dim == 1 or not given:
        raise TypeError(
            "nothing to sample: give points=(N, dim), every coordinate (x=, y=[, z=]) for "
            "points, or fix all coordinates but one for a line cut"
        )
    if dim - len(given) > 1:
        raise NotImplementedError(
            "a 3D probe with one fixed coordinate is a plane: use finest(box=...) and take "
            "the slice, or fix two coordinates for a line"
        )
    for a, v in given.items():
        if np.ndim(v) != 0:
            raise TypeError(f"{AXES[a]}={v!r}: a line cut fixes scalar coordinates; give all "
                            "coordinates as arrays for points instead")
    free = next(a for a in range(dim) if a not in given)
    # validate/wrap the fixed coordinates with a point on the line
    probe_pt = np.zeros((1, dim))
    for a, v in given.items():
        probe_pt[0, a] = float(v)
    domain_lo = 0.0 if hier.extent is None else max(hier.extent.lo[free], 0.0)
    probe_pt[0, free] = domain_lo
    probe_pt = _prepare(hier, probe_pt)
    fixed = {a: float(probe_pt[0, a]) for a in given}
    coords = line_coords(hier, fixed, free)
    pts = np.empty((len(coords), dim))
    pts[:, free] = coords
    for a, v in fixed.items():
        pts[:, a] = v
    return coords, pack(sample(snap, pts))
