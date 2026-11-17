"""Reductions over the composite AMR grid.

"Best data": every point of the extent is counted once, with the finest level available
there. Level ``l`` contributes its patches minus the pieces covered by level ``l + 1``
(``HierGeometry.covered``); the fine cells inside a covered piece and the coarse cells
outside it partition the domain exactly. The subtraction is done on node WEIGHTS (not on
sums), and coarse nodes whose cells are all covered get weight 0 and are never read into
a sum: coarse values under fine patches (even NaN) cannot leak into the result.

Integral: per cell, the multilinear (trapezoid) rule = exact integral of the multilinear
interpolant of the 2**dim corner values. Where the extent keeps only the part ``[s, e]``
(``0 <= s < e <= 1``, fractions of the cell) of a cell along an axis, the two node weights
of that cell along that axis are ``w0 = (e - s) - (e**2 - s**2) / 2`` and
``w1 = (e**2 - s**2) / 2`` times ``dl`` (exact for the linear interpolant on the kept
part); node weights are outer products over the axes. Patches and covered pieces use the
same weights. Linear fields are integrated exactly, whatever the extent.

All sums are accumulated in float64, block by block (no per-node array beyond a few
~32 MB temporaries). Nothing here reads files except through ``LevelData.buffer`` (which
reads the whole level of file-backed data first).

Two paths, same weights: patches whose cells all have full weight (inside the extent) and
that no finer patch covers use, per shape group, the outer product of cached 1-D vectors
``[h/2, h, ..., h, h/2]`` and one ``einsum`` over the group's ``[nodes, patch]`` view (no
per-node weight array); the other patches (cut by the extent, partly covered) build their
node weights block by block.

Fork safety: no BLAS call here (``@``, ``np.dot``, ``np.inner``, ``np.matmul``...). A
threaded BLAS (OpenBLAS, MKL) used before ``fork()`` hangs forever in the child process
(``multiprocessing`` "fork" start method); ``np.einsum`` without ``optimize`` uses numpy's
own loops.
"""

from __future__ import annotations

import functools
import math

import numpy as np

from ..boxes import box_flat_indices, concat_ranges
from ..boxes import snap_int as _snap_int
from .arith import SCALAR_NAMES, slot_blocks


def _bounds(hier):
    """Extent clipped to the domain: float64 ``(lo, hi)`` (dim,)."""
    domain = np.asarray(hier.domain, dtype=np.float64)
    if hier.extent is None:
        return np.zeros(hier.dim), domain
    lo, hi = hier.extent.clipped(domain)
    return np.asarray(lo, np.float64), np.asarray(hi, np.float64)


def measure(hier):
    """Volume (3D), area (2D) or length (1D) of the extent of ``hier`` inside the domain.

    Whole domain when the extent is None; 0.0 for a degenerate box (x0 == x1 on some axis).
    """
    lo, hi = _bounds(hier)
    return float(np.prod(np.clip(hi - lo, 0.0, None)))


class _AxisWeights:
    """1D trapezoid weights of one level along one axis, with the extent cut.

    Computed for the requested cells only (closed form, no table over the level's N cells:
    a deep level can have 2**27 cells per axis while its patches hold a few nodes).
    wl(c), wr(c) : weights of the left/right node of cells c (int64 array, float64 result);
                   0 for cells outside [0, N - 1] (x0, x1 lie in the domain)
    wint(i)      : weight of nodes i as interior nodes of a box: wr(i-1) + wl(i)
    cmin, cmax   : first/last cell of positive weight (cmin > cmax: none)
    fmin, fmax   : first/last cell of full weight (s == 0, e == 1; one range)
    """

    def __init__(self, ncells, dl, x0, x1):
        self.dl = dl
        self.a = float(_snap_int(x0 / dl))
        self.b = float(_snap_int(x1 / dl))
        n = int(ncells)
        # e > s  <=>  a - 1 < c < b (and a < b);  s == 0 and e == 1  <=>  a <= c <= b - 1
        if self.a < self.b:
            self.cmin = max(0, math.floor(self.a))
            self.cmax = min(n - 1, math.ceil(self.b) - 1)
        else:
            self.cmin, self.cmax = 1, 0
        self.fmin = max(0, math.ceil(self.a))
        self.fmax = min(n - 1, math.floor(self.b) - 1)
        if self.fmin > self.fmax:
            self.fmin, self.fmax = 1, 0

    def _w(self, c):
        """(wl, wr) of cells ``c`` (same float operations as a per-cell table)."""
        c = np.asarray(c, dtype=np.float64)
        s = np.clip(self.a - c, 0.0, 1.0)
        e = np.clip(self.b - c, 0.0, 1.0)
        e = np.maximum(e, s)
        w1 = (e * e - s * s) / 2
        w0 = (e - s) - w1
        return w0 * self.dl, w1 * self.dl

    def wl(self, c):
        return self._w(c)[0]

    def wr(self, c):
        return self._w(c)[1]

    def wint(self, i):
        return self.wl(i) + self.wr(i - 1)

    def box_vectors(self, glo, nnode):
        """Concatenated per-box weight vectors for boxes of nodes ``glo .. glo + nnode - 1``.

        Returns ``(v, start)``: ``v[start[b] + t]`` = weight of the t-th node of box b.
        """
        start = np.cumsum(nnode) - nnode
        v = self.wint(concat_ranges(glo, nnode))
        v[start] = self.wl(glo)
        v[start + nnode - 1] = self.wr(glo + nnode - 2)
        return v, start


#: above this many shape groups in a level, every patch takes the block path (the group
#: path is one Python iteration per group; up to 10^4-10^5 groups are possible)
GROUP_LOOP_MAX = 4096


@functools.lru_cache(maxsize=1024)
def _full_vector(n, h):
    """Trapezoid node weights of ``n`` nodes over full cells of width ``h``: h/2, h, ..., h/2.

    Equal bit for bit to ``_AxisWeights.box_vectors`` on full cells (0.5 h is exact and
    0.5 h + 0.5 h == h). Cached (keyed by (n, h)); read-only; a few hundred bytes each.
    """
    v = np.full(n, h, dtype=np.float64)
    v[0] = 0.5 * h
    v[-1] = 0.5 * h
    v.flags.writeable = False
    return v


def _group_weights(shape, dl):
    """Per-node weights of one whole full-weight patch of node ``shape``, x fastest."""
    w = _full_vector(int(shape[0]), float(dl[0]))
    for a in range(1, len(shape)):
        w = np.multiply.outer(_full_vector(int(shape[a]), float(dl[a])), w).reshape(-1)
    return w


def _group_cell_counts(shape):
    """Per-node number of the patch's cells around each node (2**dim inside), x fastest."""
    def ends(n):
        v = np.full(int(n), 2, dtype=np.int64)
        v[0] = v[-1] = 1
        return v

    c = ends(shape[0])
    for a in range(1, len(shape)):
        c = np.multiply.outer(ends(shape[a]), c).reshape(-1)
    return c


def _tiled(geom, slots, per_shape):
    """``per_shape(shape)`` repeated for every patch of ``slots`` (sorted store slots).

    Store order sorts patches by shape, so equal shapes are contiguous in ``slots``: one
    ``np.tile`` per shape run (at most one per shape group).
    """
    shape = geom.shape[slots]
    cuts = np.concatenate([[0], np.nonzero((shape[1:] != shape[:-1]).any(axis=1))[0] + 1,
                           [len(slots)]])
    return np.concatenate([np.tile(per_shape(shape[a]), b - a)
                           for a, b in zip(cuts[:-1].tolist(), cuts[1:].tolist())])


def _full_weight_patches(geom, tabs):
    """bool (npatch,): patches whose cells all have full weight (inside the extent)."""
    inside = np.ones(geom.npatch, dtype=bool)
    for a, t in enumerate(tabs):
        inside &= (geom.lower[:, a] >= t.fmin) & (geom.upper[:, a] <= t.fmax)
    return inside


def _clean_integral(geom, bufs, clean):
    """Sum over the clean patches, per component (float64), one einsum per shape group."""
    totals = np.zeros(len(bufs))
    for row in geom.groups:
        p0, p1 = int(row["p0"]), int(row["p1"])
        sel = clean[p0:p1]
        if not sel.any():
            continue
        w = _group_weights(row["shape"], geom.dl)
        n0, n1 = int(row["n0"]), int(row["n1"])
        for c, f in enumerate(bufs):
            stack = f[n0:n1].reshape((len(w), p1 - p0), order="F")
            # per-patch sums: values of the other patches of the group (NaN under a fine
            # patch, ...) stay in their own column and are dropped by ``sel``
            per_patch = np.einsum("nk,n->k", stack, w)
            totals[c] += per_patch[sel].sum()
    return totals


def _weights_tables(geom, hier):
    lo, hi = _bounds(hier)
    return [_AxisWeights(int(geom.domain_cells[a]), float(geom.dl[a]), lo[a], hi[a])
            for a in range(geom.dim)]


def _crop(tabs, clo, chi):
    """Crop global cell boxes to the cells of positive weight; returns (clo, chi, keep)."""
    cmin = np.array([t.cmin for t in tabs], dtype=np.int64)
    cmax = np.array([t.cmax for t in tabs], dtype=np.int64)
    clo = np.maximum(clo, cmin)
    chi = np.minimum(chi, cmax)
    return clo, chi, (clo <= chi).all(axis=1)


def _ends_vectors(nnode, inside, end, dtype):
    """Per-box vectors equal to ``inside`` except ``end`` at both ends (concatenated)."""
    start = np.cumsum(nnode) - nnode
    v = np.full(int(nnode.sum()), inside, dtype=dtype)
    v[start] = end
    v[start + nnode - 1] = end
    return v, start


def _outer(vectors, nnode):
    """Per-node outer product of per-box axis vectors, nodes of every box x fastest.

    vectors : one ``(v, start)`` per axis: ``v[start[b] + t]`` is the factor of the t-th node
              of box b along that axis. nnode : (n, dim) node counts of the boxes.
    Returns the products for the boxes one after the other (the order of
    ``box_flat_indices`` with run_axis=0).
    """
    nb, dim = nnode.shape
    m = np.prod(nnode[:, 1:], axis=1) if dim > 1 else np.ones(nb, np.int64)
    line_box = np.repeat(np.arange(nb, dtype=np.int64), m)
    line_len = nnode[line_box, 0]
    vx, sx = vectors[0]
    out = vx[concat_ranges(sx[line_box], line_len)]
    if dim > 1:
        local = concat_ranges(np.zeros(nb, np.int64), m)
        ny = nnode[line_box, 1]
        vy, sy = vectors[1]
        wline = vy[sy[line_box] + local % ny]
        if dim > 2:
            vz, sz = vectors[2]
            wline = wline * vz[sz[line_box] + local // ny]
        out *= np.repeat(wline, line_len)
    return out


def _box_weights(tabs, glo, nnode):
    """Integration weight of every node of boxes of whole cells (cut by the extent)."""
    return _outer([t.box_vectors(glo[:, a], nnode[:, a]) for a, t in enumerate(tabs)], nnode)


def _box_cell_counts(nnode):
    """Number of the box's cells adjacent to every node (2**dim inside, fewer on faces)."""
    vec = [_ends_vectors(nnode[:, a], 2, 1, np.int64) for a in range(nnode.shape[1])]
    return _outer(vec, nnode)


# the level sum is not sum_patches T(f) - sum_pieces T(f): the piece WEIGHTS are
# subtracted node by node instead (same value in exact arithmetic), so that
# NaN or garbage under fine patches cannot leak in through cancellation.
def _level_integral(ld, hier):
    """Contribution of one level, one float64 value per component.

    Patches are cropped to the cells of positive weight (extent). Per block of patches,
    the node weights are the patch weights minus the weights of the covered pieces in
    them (the cells of a piece belong to the finer level). Nodes all of whose cells are
    covered (counted exactly with integers) get weight 0 and their values are not read
    into the sum, so whatever the coarse level holds under a fine patch (coarsened data,
    garbage, NaN) never contributes.
    """
    geom = ld.geom
    tabs = _weights_tables(geom, hier)
    bufs = list(ld.buffers().values())
    cover = hier.covered(geom.level)
    totals = np.zeros(len(bufs))
    clean = np.zeros(geom.npatch, dtype=bool)
    by_group = len(geom.groups) <= GROUP_LOOP_MAX
    inside = None  # full-weight patches, only computed (and used) when by_group
    if by_group:
        inside = _full_weight_patches(geom, tabs)
        clean = inside.copy()
        clean[cover.slot] = False
        totals += _clean_integral(geom, bufs, clean)
    clo, chi, keep = _crop(tabs, geom.lower, geom.upper)
    keep &= ~clean
    slots = np.nonzero(keep)[0]
    clo, chi = clo[keep], chi[keep]
    nnode = chi - clo + 2
    whole = len(slots) == geom.npatch and np.array_equal(clo, geom.lower) and np.array_equal(
        chi, geom.upper)
    plo = cover.lo + geom.lower[cover.slot]
    phi = cover.hi + geom.lower[cover.slot]
    plo, phi, pkeep = _crop(tabs, plo, phi)
    pbox = np.searchsorted(slots, cover.slot[pkeep])  # index of the piece's patch in ``slots``
    plo, phi = plo[pkeep], phi[pkeep]
    for b0, b1 in slot_blocks(np.prod(nnode, axis=1)):
        s, g, n = slots[b0:b1], clo[b0:b1], nnode[b0:b1]
        # whole full-weight patches (here: covered ones) tile per-shape weights
        uncut = by_group and bool(inside[s].all())
        if uncut:
            w = _tiled(geom, s, lambda sh: _group_weights(sh, geom.dl))
        else:
            w = _box_weights(tabs, g, n)
        full = None
        k0, k1 = np.searchsorted(pbox, [b0, b1])
        if k1 > k0:
            # piece nodes as positions in this block's node enumeration
            start = np.cumsum(np.prod(n, axis=1)) - np.prod(n, axis=1)
            k = pbox[k0:k1] - b0
            pn = phi[k0:k1] - plo[k0:k1] + 2
            pos = box_flat_indices(start[k], n[k], plo[k0:k1] - g[k], pn)
            w -= np.bincount(pos, _box_weights(tabs, plo[k0:k1], pn), minlength=len(w))
            covered_cells = np.bincount(pos, _box_cell_counts(pn), minlength=len(w))
            cells = _tiled(geom, s, _group_cell_counts) if uncut else _box_cell_counts(n)
            full = covered_cells == cells
            w[full] = 0.0
        if whole:
            idx = slice(int(geom.offset[s[0]]), int(geom.offset[s[-1] + 1]))
        elif uncut:
            idx = concat_ranges(geom.offset[s], geom.size[s])
        else:
            idx = box_flat_indices(geom.offset[s], geom.shape[s], g - geom.lower[s], n)
        for c, f in enumerate(bufs):
            vals = f[idx].astype(np.float64)
            if full is not None:
                vals[full] = 0.0
            totals[c] += np.einsum("i,i->", vals, w)  # not ``@``: no BLAS (fork safety)
    return totals


def _require_snapshot(snap, what):
    from ..store import Snapshot

    if not isinstance(snap, Snapshot):
        raise TypeError(f"{what} expects an AMR field snapshot, got {type(snap).__name__}")


def _shape_result(values, names):
    values = np.asarray(values, dtype=np.float64)
    return float(values[0]) if tuple(names) == SCALAR_NAMES else values


def integral(snap):
    """Integral of the field over its extent, with the best data at every point.

    Returns a float for a scalar field, a float64 array ``(ncomp,)`` (components in
    ``snap.names`` order) otherwise. Units: field units times length**dim. Exact for fields
    linear in each coordinate. Reads every level fully (file-backed data).
    """
    _require_snapshot(snap, "integral")
    total = np.zeros(len(snap.names))
    for ld in snap.levels:
        total += _level_integral(ld, snap.hier)
    return _shape_result(total, snap.names)


def mean(snap):
    """``integral(snap) / measure(snap.hier)`` (area-weighted mean, best data).

    Raises ValueError when the extent has zero measure.
    """
    _require_snapshot(snap, "mean")
    vol = measure(snap.hier)
    if vol <= 0:
        raise ValueError("mean: the field's box has zero length/area/volume; widen the box")
    total = integral(snap)
    return total / vol


def _node_range(geom, hier):
    """Global node index range ``[nlo, nhi]`` (dim,) inside the CLOSED extent."""
    lo, hi = _bounds(hier)
    dl = np.asarray(geom.dl)
    nlo = np.ceil(_snap_int(lo / dl)).astype(np.int64)
    nhi = np.floor(_snap_int(hi / dl)).astype(np.int64)
    return nlo, nhi


# a box lying between two node planes holds no node: min/max raise ValueError
# (the integral and mean of such a box are still defined).
def _extreme(snap, how, finite=False):
    _require_snapshot(snap, how)
    # fmin/fmax ignore NaN: a field masked at some nodes (e.g. divB at coarse-fine faces)
    # still has a min and a max; the result is NaN only if every kept node is NaN
    reduce = np.fmin if how == "min" else np.fmax
    best = [None] * len(snap.names)
    has_extent = snap.hier.extent is not None
    for ld in snap.levels:
        geom = ld.geom
        nlo, nhi = _node_range(geom, snap.hier)
        slots = None
        if has_extent:  # read only the patches with a node in the extent (_block_mask's rule)
            inside = (np.maximum(geom.lower, nlo) <= np.minimum(geom.upper + 1, nhi)).all(1)
            slots = np.flatnonzero(inside)
        loaded = ld.ensure(slots)
        if loaded is None:  # a level without patches
            continue
        bufs = [loaded[n] for n in ld.names]
        cover = _covered_node_boxes(geom, snap.hier.covered(geom.level))
        for i0, i1 in slot_blocks(geom.size):
            n0, n1 = int(geom.offset[i0]), int(geom.offset[i1])
            keep = _block_mask(geom, has_extent, nlo, nhi, cover, i0, i1)
            if keep is not None and not keep.any():
                continue
            for k, f in enumerate(bufs):
                vals = f[n0:n1] if keep is None else f[n0:n1][keep]
                if finite:
                    vals = vals[np.isfinite(vals)]
                    if not vals.size:
                        continue
                v = reduce.reduce(vals)
                best[k] = v if best[k] is None else reduce(best[k], v)
    if all(v is None for v in best):
        raise ValueError(f"{how}: no node of the field lies inside its box")
    # None: every kept node of that component is non-finite (finite=True)
    return _shape_result([np.nan if v is None else float(v) for v in best], snap.names)


def _covered_node_boxes(geom, cover):
    """Closed covered node boxes cut for EVERY patch holding a copy of those nodes.

    A node on the edge of a covered region is stored by each patch that touches it, not only
    by the patch owning the covered cells; the finer level holds it, so every copy goes.
    Returns ``(slot, lo, size)`` sorted by slot: ``lo`` local to the patch, ``size`` in nodes.
    """
    dim = geom.dim
    if not len(cover):
        z = np.zeros(0, dtype=np.int64)
        return z, np.zeros((0, dim), dtype=np.int64), np.zeros((0, dim), dtype=np.int64)
    glo = cover.lo + geom.lower[cover.slot]
    ghi = cover.hi + geom.lower[cover.slot]
    pairs = geom.bin_index().pairs(glo - 1, ghi + 1)  # cells whose node box touches a piece
    a, b = pairs.a, pairs.b
    lo = np.maximum(glo[a], geom.lower[b])
    hi = np.minimum(ghi[a] + 1, geom.upper[b] + 1)
    ok = (lo <= hi).all(axis=1)
    a, b, lo, hi = a[ok], b[ok], lo[ok], hi[ok]
    order = np.argsort(b, kind="stable")
    b, lo, hi = b[order], lo[order], hi[order]
    return b, lo - geom.lower[b], hi - lo + 1


def _block_mask(geom, has_extent, nlo, nhi, cover, i0, i1):
    """Bool mask of the nodes of patches i0..i1-1 to keep, or None for all of them."""
    cslot, clo, csize = cover
    n0, n1 = int(geom.offset[i0]), int(geom.offset[i1])
    p0, p1 = np.searchsorted(cslot, [i0, i1])
    if not has_extent and p0 == p1:
        return None
    keep = np.ones(n1 - n0, dtype=bool)
    if has_extent:
        keep[:] = False
        lo = np.maximum(geom.lower[i0:i1], nlo)
        hi = np.minimum(geom.upper[i0:i1] + 1, nhi)
        inside = (lo <= hi).all(axis=1)
        s = np.arange(i0, i1)[inside]
        idx = box_flat_indices(geom.offset[s], geom.shape[s], (lo - geom.lower[i0:i1])[inside],
                               (hi - lo + 1)[inside])
        keep[idx - n0] = True
    if p1 > p0:
        s = cslot[p0:p1]
        idx = box_flat_indices(geom.offset[s], geom.shape[s], clo[p0:p1], csize[p0:p1])
        keep[idx - n0] = False
    return keep


def minimum(snap, finite=False):
    """Smallest value over the best data inside the (closed) extent.

    Nodes of a level inside or ON THE BOUNDARY of a piece covered by the next level are
    skipped (the finer level holds those points). NaN values are ignored (the result is NaN
    only when every kept node is NaN); ``finite=True`` ignores +-inf as well. Float for a
    scalar field, float64 ``(ncomp,)`` otherwise. Raises ValueError if no node is left.
    """
    return _extreme(snap, "min", finite)


def maximum(snap, finite=False):
    """Largest value over the best data inside the (closed) extent (see ``minimum``)."""
    return _extreme(snap, "max", finite)


# per_level ignores the extent beyond the whole patches present (master's meaning).
def per_level(snap, how):
    """``{level number: value}`` with master's per-level meaning (no weights, no coverage).

    how : "mean", "min", "max" or "sum" over all nodes of each level's own data (nodes shared
          by two patches count twice; the extent selects whole patches only). float64; a
          float per level for a scalar field, ``(ncomp,)`` arrays otherwise.
    min/max ignore NaN as ``minimum``/``maximum`` do (NaN only when every node of the
    level is NaN); mean/sum propagate it.
    """
    _require_snapshot(snap, "per_level")
    funcs = {
        "mean": lambda f: np.mean(f, dtype=np.float64),
        "sum": lambda f: np.sum(f, dtype=np.float64),
        "min": np.fmin.reduce,
        "max": np.fmax.reduce,
    }
    if how not in funcs:
        raise ValueError(f"per_level: how must be one of {tuple(funcs)}, got {how!r}")
    return {ld.level: _shape_result([float(funcs[how](b)) for b in ld.buffers().values()],
                                    snap.names)
            for ld in snap.levels}
