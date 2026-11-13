"""Level-aware finite differences on node-centred AMR data (DESIGN 2 "derivatives", 7.2).

Every derivative uses data of ONE level only (never coarse data):
  * interior nodes of a patch: centred differences, done as one flat shift per stride run
    of the level (``LevelGeometry.runs``), not per patch;
  * the first and last node layers of every patch along the axis: one-sided formulas
    (3-point, second order, for the first derivative; 2 nodes along the axis: 2-point,
    first order);
  * then, for every face node that has a same-level neighbour holding the node one step
    outside (another patch touching it, or a periodic image), the centred formula is
    written over the one-sided value (``LevelGeometry.faces`` -> ``FaceTable``).
A node shared by two patches therefore receives the same value from both when both use
the centred formula (at a one-sided end, two patches of different widths along the axis,
e.g. 2 nodes and 5 nodes, use different one-sided formulas and may differ). One-sided values remain at domain boundaries (non periodic) and at
coarse-fine faces: there the derivative is less accurate, a "noise seam" visible on noisy
data (``one_sided_mask`` / ``div(mask_faces=True)`` show where).

Arrays are flat, in store order, x fastest inside a patch; the dtype of the field is kept
for floating types (float32 stays float32; integers give float64). ``h`` is the mesh size
of the level along the axis.
"""

from __future__ import annotations

import numpy as np

from ..boxes import box_flat_indices, fortran_strides
from ..plan import ReadPlan
from ..store import LevelData, Snapshot
from . import arith
from .arith import SCALAR_NAMES, VECTOR_NAMES, slot_blocks


def _float_dtype(dtype):
    dtype = np.dtype(dtype)
    return dtype if np.issubdtype(dtype, np.floating) else np.dtype(np.float64)


def _field(ld, name):
    f = ld.buffer(name)
    dt = _float_dtype(f.dtype)
    return (f if f.dtype == dt else f.astype(dt)), dt


def _check_axis(axis):
    if axis not in (0, 1, 2):
        raise ValueError(f"axis must be 0, 1 or 2 (x, y, z), got {axis!r}")


def end_layer_blocks(geom, axis, max_nodes=None):
    """First/last node layers along ``axis``, by blocks of consecutive patches.

    Yields ``(lo, hi, stride, count)`` int64 arrays (one entry per node of the layers of the
    block's patches, like ``LevelGeometry.end_layers`` restricted to the block).
    """
    shape = geom.shape
    size = shape.copy()
    size[:, axis] = 1
    per_patch = np.prod(size, axis=1)
    strides = fortran_strides(shape)[:, axis]
    for p0, p1 in slot_blocks(per_patch, max_nodes):
        lo = box_flat_indices(geom.offset[p0:p1], shape[p0:p1],
                              np.zeros((p1 - p0, geom.dim), np.int64), size[p0:p1])
        stride = np.repeat(strides[p0:p1], per_patch[p0:p1])
        count = np.repeat(shape[p0:p1, axis], per_patch[p0:p1])
        yield lo, lo + (count - 1) * stride, stride, count


def _face_blocks(table):
    max_rows = arith.BLOCK_NODES
    for c in range(0, len(table), max_rows):
        sl = slice(c, c + max_rows)
        sign = table.sign[sl]
        inner = table.inner_rows(c, min(c + max_rows, len(table)))  # computed, once per block
        plus = np.where(sign > 0, table.outer[sl], inner)
        minus = np.where(sign > 0, inner, table.outer[sl])
        yield table.face[sl], plus, minus


# --------------------------------------------------------------------------- level kernels


# DESIGN-NOTE: an axis the data does not have returns zeros (grad/curl need d/dz = 0 in 2D)
# instead of raising.
def ddx(ld, name, axis):
    """First derivative of component ``name`` of LevelData ``ld`` along ``axis`` (0, 1, 2).

    Returns a new flat array ``(nnodes,)`` in store order (x fastest per patch), in the
    field's floating dtype. With ``h = dl[axis]`` of the level:
      interior        (f[i+1] - f[i-1]) / (2h)
      first layer     (-3 f[0] + 4 f[1] - f[2]) / (2h)      (3+ nodes along the axis)
      last layer      ( 3 f[n] - 4 f[n-1] + f[n-2]) / (2h)
      2-node patches  (f[1] - f[0]) / h at both layers (first order)
      face node with a same-level neighbour (FaceTable): (f[+h] - f[-h]) / (2h), centred.
    An axis the data does not have (e.g. z in 2D) gives zeros. Reads the whole level first
    (file-backed data). Raises ValueError for a bad axis, KeyError for a bad name.
    """
    _check_axis(axis)
    geom = ld.geom
    f, dt = _field(ld, name)
    if axis >= geom.dim:
        return np.zeros(geom.nnodes, dtype=dt)
    d = np.empty(geom.nnodes, dtype=dt)
    two_h = dt.type(2.0 * geom.dl[axis])
    h = dt.type(geom.dl[axis])
    # interior: one flat shift per stride run (end layers get wrong values, fixed below)
    for s, a, b in zip(*geom.runs(axis)):
        s, a, b = int(s), int(a), int(b)
        if b - a > 2 * s:
            out = d[a + s : b - s]
            np.subtract(f[a + 2 * s : b], f[a : b - 2 * s], out=out)
            np.divide(out, two_h, out=out)
    # end layers: one-sided
    for lo, hi, s, n in end_layer_blocks(geom, axis):
        big = n >= 3
        lo3, hi3, s3 = lo[big], hi[big], s[big]
        d[lo3] = (-3 * f[lo3] + 4 * f[lo3 + s3] - f[lo3 + 2 * s3]) / two_h
        d[hi3] = (3 * f[hi3] - 4 * f[hi3 - s3] + f[hi3 - 2 * s3]) / two_h
        lo2, hi2 = lo[~big], hi[~big]
        d2 = (f[hi2] - f[lo2]) / h
        d[lo2] = d2
        d[hi2] = d2
    # faces with a same-level neighbour: centred again
    for face, plus, minus in _face_blocks(geom.faces(axis)):
        d[face] = (f[plus] - f[minus]) / two_h
    return d


# DESIGN-NOTE: DESIGN 7.2 says "< 4 nodes: 3-point" without covering 2 nodes; a 2-node patch
# end without a neighbour has no curvature information and gets 0.
def d2dx2(ld, name, axis):
    """Second derivative of component ``name`` of ``ld`` along ``axis`` (flat, store order).

    With ``h = dl[axis]``:
      interior          (f[i+1] - 2 f[i] + f[i-1]) / h**2
      end layers (4+ nodes along the axis): (2 f[0] - 5 f[1] + 4 f[2] - f[3]) / h**2
                        (second order, exact for cubics), mirrored at the last layer
      3 nodes           both end layers take the value of the middle node (3-point)
      2 nodes           0 (no curvature information)
      face node with a same-level neighbour: (f[+h] - 2 f[face] + f[-h]) / h**2
    Missing axis -> zeros. Same reading and errors as ``ddx``.
    """
    _check_axis(axis)
    geom = ld.geom
    f, dt = _field(ld, name)
    if axis >= geom.dim:
        return np.zeros(geom.nnodes, dtype=dt)
    d = np.empty(geom.nnodes, dtype=dt)
    h2 = dt.type(geom.dl[axis] * geom.dl[axis])
    for s, a, b in zip(*geom.runs(axis)):
        s, a, b = int(s), int(a), int(b)
        # chunks bound the temporary of 2 f
        step = arith.BLOCK_NODES
        for q0 in range(a + s, b - s, step):
            q1 = min(q0 + step, b - s)
            out = d[q0:q1]
            tmp = 2 * f[q0:q1]
            np.subtract(f[q0 + s : q1 + s], tmp, out=out)
            np.add(out, f[q0 - s : q1 - s], out=out)
            np.divide(out, h2, out=out)
    for lo, hi, s, n in end_layer_blocks(geom, axis):
        four = n >= 4
        l4, h4, s4 = lo[four], hi[four], s[four]
        d[l4] = (2 * f[l4] - 5 * f[l4 + s4] + 4 * f[l4 + 2 * s4] - f[l4 + 3 * s4]) / h2
        d[h4] = (2 * f[h4] - 5 * f[h4 - s4] + 4 * f[h4 - 2 * s4] - f[h4 - 3 * s4]) / h2
        three = n == 3
        l3, s3 = lo[three], s[three]
        mid = (f[l3 + 2 * s3] - 2 * f[l3 + s3] + f[l3]) / h2
        d[l3] = mid
        d[hi[three]] = mid
        d[lo[n == 2]] = 0
        d[hi[n == 2]] = 0
    for face, plus, minus in _face_blocks(geom.faces(axis)):
        d[face] = (f[plus] - 2 * f[face] + f[minus]) / h2
    return d


def one_sided_mask(geom):
    """Bool flat array ``(nnodes,)``: nodes where at least one axis used a one-sided formula.

    That is: nodes on the first/last layer of their patch along some axis, without a
    same-level neighbour across that face (domain boundary or coarse-fine face).
    """
    mask = np.zeros(geom.nnodes, dtype=bool)
    for axis in range(geom.dim):
        m = np.zeros(geom.nnodes, dtype=bool)
        for lo, hi, _, _ in end_layer_blocks(geom, axis):
            m[lo] = True
            m[hi] = True
        m[geom.faces(axis).face] = False
        mask |= m
    return mask


# --------------------------------------------------------------------------- snapshot operators


def _require(snap, names, what, expected):
    if not isinstance(snap, Snapshot) or tuple(snap.names) != names:
        got = snap.names if isinstance(snap, Snapshot) else type(snap).__name__
        raise TypeError(f"{what} expects {expected} (components {names}), got {got}")


def _check_result(snap, ncomp, what):
    plan = ReadPlan()
    for ld in snap.levels:
        plan.add_read(ld)
    dt = _float_dtype(snap.levels[0].dtype)
    plan.add_result(sum(g.nnodes for g in snap.hier.levels), ncomp, dt, what=what)
    plan.check(what=what)


def _new(snap, per_level):
    levels = tuple(LevelData.from_arrays(ld.geom, ld.time, arrays)
                   for ld, arrays in zip(snap.levels, per_level))
    return Snapshot(snap.hier, levels)


def grad(snap):
    """Gradient of a scalar Snapshot -> vector Snapshot ("x", "y", "z").

    Components along axes the data does not have (y, z in 1D; z in 2D) are zeros.
    Raises TypeError for a non-scalar input.
    """
    _require(snap, SCALAR_NAMES, "grad", "a scalar field")
    _check_result(snap, 3, "grad")
    return _new(snap, [{c: ddx(ld, "value", a) for a, c in enumerate(VECTOR_NAMES)}
                       for ld in snap.levels])


def _diff(p, m, n, dt):
    """p - m where None means 0."""
    if p is None and m is None:
        return np.zeros(n, dtype=dt)
    if m is None:
        return p
    if p is None:
        np.negative(m, out=m)
        return m
    np.subtract(p, m, out=p)
    return p


def curl(snap):
    """Curl of a vector Snapshot -> vector Snapshot, with the derivatives of ``ddx``.

    Derivatives along missing axes are 0: a 2D field gives all 3 components
    (``curl_z = dBy/dx - dBx/dy``, ``curl_x = dBz/dy``, ``curl_y = -dBz/dx``).
    On node data this equals PHARE's Yee J with an extra (1, 2, 1)/4 smoothing along each
    derivative direction (centred differences of node-averaged B): on under-resolved
    current sheets its peak is lower than the code's J (7-50 % lower on the sample runs).
    The Yee J itself needs Yee-centred input (phareh5 files, ``Run.GetJ(stencil="yee")``).
    Raises TypeError for a non-vector input.
    """
    _require(snap, VECTOR_NAMES, "curl", "a vector field")
    _check_result(snap, 3, "curl")
    out = []
    for ld in snap.levels:
        dim, n = ld.geom.dim, ld.geom.nnodes
        dt = _float_dtype(ld.dtype)
        ld.ensure()  # the three components in one read

        def d(name, axis, ld=ld, dim=dim):
            return ddx(ld, name, axis) if axis < dim else None

        out.append({
            "x": _diff(d("z", 1), d("y", 2), n, dt),
            "y": _diff(d("x", 2), d("z", 0), n, dt),
            "z": _diff(d("y", 0), d("x", 1), n, dt),
        })
    return _new(snap, out)


def div(snap, mask_faces=False):
    """Divergence of a vector Snapshot -> scalar Snapshot: ``sum_a ddx(B_a, a)``.

    On PHARE node data this is exactly the mean of PHARE's Yee divB over the 2**dim cells
    around each node, at nodes where every derivative is centred (away from coarse-fine
    faces and non-periodic boundaries). Where a one-sided formula was used the value is
    dominated by its truncation error on noisy data (tens to a hundred times the Yee divB
    on the sample runs); ``mask_faces=True`` sets those nodes to NaN (``min()`` and
    ``max()`` of the result ignore them; ``mean()`` and ``integral()`` return NaN). Raises TypeError for a non-vector
    input.
    """
    _require(snap, VECTOR_NAMES, "div", "a vector field")
    _check_result(snap, 1, "div")
    out = []
    for ld in snap.levels:
        ld.ensure(None, VECTOR_NAMES[: ld.geom.dim])  # one read for the components used
        r = ddx(ld, "x", 0)
        for a in range(1, ld.geom.dim):
            r += ddx(ld, VECTOR_NAMES[a], a)
        if mask_faces:
            r[one_sided_mask(ld.geom)] = np.nan
        out.append({"value": r})
    return _new(snap, out)


def laplacian(snap):
    """Laplacian ``sum_a d2dx2(f, a)``: scalar -> scalar, vector -> vector (per component).

    Raises TypeError for a Snapshot that is neither scalar nor vector.
    """
    if not isinstance(snap, Snapshot) or tuple(snap.names) not in (SCALAR_NAMES, VECTOR_NAMES):
        got = snap.names if isinstance(snap, Snapshot) else type(snap).__name__
        raise TypeError(f"laplacian expects a scalar or vector field, got {got}")
    _check_result(snap, len(snap.names), "laplacian")
    out = []
    for ld in snap.levels:
        ld.ensure()  # every component in one read
        comps = {}
        for name in ld.names:
            r = d2dx2(ld, name, 0)
            for a in range(1, ld.geom.dim):
                r += d2dx2(ld, name, a)
            comps[name] = r
        out.append(comps)
    return _new(snap, out)
