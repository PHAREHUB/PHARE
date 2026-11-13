"""Integer box tables and the vectorised index arithmetic built on them.

Conventions used everywhere in this module
------------------------------------------
* A *box table* is a pair ``(lower, upper)`` of int64 arrays of shape ``(n, dim)``:
  box ``i`` is the inclusive integer box ``lower[i] .. upper[i]`` (both ends included).
  Most of the time the integers are cell indices of one AMR level; the functions here do
  not care whether they are cells or nodes, they only see integers.
* A patch with cell box ``[lower, upper]`` has nodes ``lower .. upper + 1``: node shape
  ``upper - lower + 2`` per axis. Its values are stored x fastest (Fortran order) in a
  flat buffer starting at a node *offset*: node ``(i, j, k)`` (local, 0-based) is at
  ``offset + i + nx * j + nx * ny * k``.
* Every index array returned is int64.

Nothing here reads files; nothing loops in Python over boxes (only over bounded blocks).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import NamedTuple

import numpy as np

AXES = ("x", "y", "z")


# --------------------------------------------------------------------------- basics


def as_box_table(lower, upper):
    """Return ``(lower, upper)`` as int64 arrays of shape ``(n, dim)``.

    A single box given as two ``(dim,)`` sequences becomes a table with one row.
    Raises ValueError when the shapes differ or a box is empty (``upper < lower``).
    """
    lower = np.atleast_2d(np.asarray(lower, dtype=np.int64))
    upper = np.atleast_2d(np.asarray(upper, dtype=np.int64))
    if lower.shape != upper.shape or lower.ndim != 2:
        raise ValueError(f"box table: lower {lower.shape} and upper {upper.shape} must be (n, dim)")
    if lower.size and (upper < lower).any():
        bad = int(np.nonzero((upper < lower).any(axis=1))[0][0])
        raise ValueError(f"box {bad} is empty: lower {lower[bad]} > upper {upper[bad]}")
    return lower, upper


def coarsen(lower, upper, ratio=2):
    """Coarse cell boxes holding the fine cell boxes ``[lower, upper]`` (floor division).

    Example (ratio 2): fine cells ``[4, 9]`` -> coarse cells ``[2, 4]``;
    fine ``[5, 8]`` -> ``[2, 4]``
    (a partly covered coarse cell is included). See ``coarsen_inner`` for the opposite rule.
    """
    return np.floor_divide(lower, ratio), np.floor_divide(upper, ratio)


def coarsen_inner(lower, upper, ratio=2):
    """Coarse cell boxes *fully* covered by the fine cell boxes ``[lower, upper]``.

    Equal to ``coarsen`` for ratio-aligned boxes (lower even, upper odd for ratio 2).
    Example (ratio 2): fine ``[4, 9]`` -> ``[2, 4]``; fine ``[5, 8]`` -> ``[3, 3]``; fine ``[5, 5]``
    -> empty (returned as lower > upper, callers drop such rows).
    """
    lo = -np.floor_divide(-np.asarray(lower), ratio)  # ceil(lower / ratio)
    hi = np.floor_divide(np.asarray(upper) + 1, ratio) - 1
    return lo, hi


def refine(lower, upper, ratio=2):
    """Fine cell boxes covering the coarse cell boxes ``[lower, upper]``.

    Example (ratio 2): coarse ``[2, 4]`` -> fine ``[4, 9]``.
    """
    return np.asarray(lower) * ratio, (np.asarray(upper) + 1) * ratio - 1


def intersect(alo, ahi, blo, bhi):
    """Row-wise intersection of two box tables of equal length.

    Returns ``(lo, hi, nonempty)``: ``nonempty[i]`` is True when boxes ``i`` share at least
    one integer point (closed boxes).
    """
    lo = np.maximum(alo, blo)
    hi = np.minimum(ahi, bhi)
    return lo, hi, (lo <= hi).all(axis=-1)


def box_counts(lo, hi):
    """Number of integer points of every box, int64 ``(n,)``; 0 for empty boxes."""
    ext = np.asarray(hi, dtype=np.int64) - np.asarray(lo, dtype=np.int64) + 1
    return np.prod(np.clip(ext, 0, None), axis=-1)


def periodic_shifts(domain_cells, periodic, include_zero=True):
    """All shift vectors ``(k, dim)`` int64 made of ``{-N_a, 0, +N_a}`` on the periodic axes.

    ``domain_cells`` (dim,) is the number of cells of the level along each axis. The zero
    shift comes first when ``include_zero``. Example: 2D, periodic in x only, N = (8, 4):
    ``[[0, 0], [-8, 0], [8, 0]]``.
    """
    domain_cells = np.asarray(domain_cells, dtype=np.int64)
    periodic = np.asarray(periodic, dtype=bool)
    choices = [(0, -int(n), int(n)) if p else (0,) for n, p in zip(domain_cells, periodic)]
    grids = np.meshgrid(*[np.array(c, dtype=np.int64) for c in choices], indexing="ij")
    shifts = np.stack([g.ravel() for g in grids], axis=1)
    if not include_zero:
        shifts = shifts[(shifts != 0).any(axis=1)]
    return shifts


# --------------------------------------------------------------------------- flat indices


def concat_ranges(starts, lengths):
    """Concatenation of ``arange(s, s + n)`` for every ``(s, n)``, without a Python loop.

    starts, lengths : int arrays of equal length; lengths >= 0.
    Returns int64 ``(sum(lengths),)``. Example: ``concat_ranges([10, 0], [3, 2])`` ->
    ``[10, 11, 12, 0, 1]``.
    """
    starts = np.asarray(starts, dtype=np.int64)
    lengths = np.asarray(lengths, dtype=np.int64)
    total = int(lengths.sum())
    if total == 0:
        return np.zeros(0, dtype=np.int64)
    # each element = its range start + (its position in the output - where its range begins)
    begins = np.cumsum(lengths) - lengths
    return np.repeat(starts - begins, lengths) + np.arange(total, dtype=np.int64)


def fortran_strides(shapes):
    """Element strides ``(n, dim)`` int64 of x-fastest arrays of the given node shapes."""
    shapes = np.asarray(shapes, dtype=np.int64)
    strides = np.ones_like(shapes)
    if shapes.shape[1] > 1:
        strides[:, 1:] = np.cumprod(shapes[:, :-1], axis=1)
    return strides


def box_flat_indices(offsets, shapes, lo, size, run_axis=0):
    """Flat buffer indices of the nodes of sub-boxes of patches, x fastest.

    Each row ``i`` describes one sub-box of one patch:
      offsets : (n,) int64   node offset of the patch in the level's flat buffer
      shapes  : (n, dim)     node shape of that patch (``upper - lower + 2``)
      lo      : (n, dim)     first node of the sub-box, LOCAL to the patch (0-based)
      size    : (n, dim)     number of nodes of the sub-box along each axis (>= 0)
    The sub-box must lie inside its patch (not checked; out-of-range input gives wrong indices).

    Returns int64 ``(sum(prod(size, axis=1)),)``: the sub-boxes one after the other in row
    order, each one x fastest (Fortran order). For one row this is exactly::

        offsets[i] + np.arange(prod(shapes[i])).reshape(shapes[i], order="F")[
            lo[i,0]:lo[i,0]+size[i,0], lo[i,1]:lo[i,1]+size[i,1]].ravel(order="F")

    Example (2D): a patch of node shape (4, 3) at offset 100, sub-box lo=(1, 1), size=(2, 2)
    -> ``[105, 106, 109, 110]`` (node (1,1) is 100 + 1 + 4*1).
    Cost O(number of output indices); use ``iter_box_flat_indices`` to bound memory.

    run_axis : axis enumerated fastest (default 0 = the x-fastest order above). Another axis,
               or "auto" (the axis with the most nodes in total), gives the SAME set of indices
               per sub-box in another order (that axis fastest, then the others in increasing
               order); use it when the order does not matter: it is several times faster
               on sub-boxes one node thick along x (faces along x). The order depends only on
               ``size``, so two calls with equal ``size`` enumerate nodes in matching order.
    """
    offsets = np.asarray(offsets, dtype=np.int64)
    shapes = np.atleast_2d(np.asarray(shapes, dtype=np.int64))
    lo = np.atleast_2d(np.asarray(lo, dtype=np.int64))
    size = np.clip(np.atleast_2d(np.asarray(size, dtype=np.int64)), 0, None)
    n, dim = shapes.shape
    if n == 0:
        return np.zeros(0, dtype=np.int64)
    strides = fortran_strides(shapes)
    if run_axis == "auto":
        run_axis = int(np.argmax(size.sum(axis=0)))
    others = [a for a in range(dim) if a != run_axis]
    # enumerate rows (all axes but the run axis), then expand every row into a strided run
    nrows = np.prod(size[:, others], axis=1) if others else np.ones(n, dtype=np.int64)
    nrows = np.where(size[:, run_axis] > 0, nrows, 0)
    base = offsets + (lo * strides).sum(axis=1)
    piece = np.repeat(np.arange(n, dtype=np.int64), nrows)
    row = concat_ranges(np.zeros(n, dtype=np.int64), nrows)
    start = base[piece]
    for a in others:
        s = size[piece, a]
        start += (row % s) * strides[piece, a]
        row //= s
    length = size[piece, run_axis]
    if run_axis == 0:  # contiguous runs
        return concat_ranges(start, length)
    step = np.repeat(strides[piece, run_axis], length)
    k = concat_ranges(np.zeros(len(start), dtype=np.int64), length)
    return np.repeat(start, length) + k * step


def iter_box_flat_indices(offsets, shapes, lo, size, max_nodes=16 * 2**20):
    """``box_flat_indices`` in blocks of whole sub-boxes holding at most ``max_nodes`` nodes.

    Yields ``(i0, i1, idx)``: rows ``i0 .. i1-1`` and their concatenated flat indices, so that
    concatenating every ``idx`` gives ``box_flat_indices(...)``. A single sub-box larger
    than ``max_nodes`` is yielded alone. Loops once per block, never per row.
    """
    size = np.clip(np.atleast_2d(np.asarray(size, dtype=np.int64)), 0, None)
    counts = np.prod(size, axis=1)
    ends = np.cumsum(counts)
    n = len(counts)
    i0 = 0
    while i0 < n:
        before = ends[i0] - counts[i0]
        i1 = int(np.searchsorted(ends, before + max_nodes, side="right"))
        i1 = max(i1, i0 + 1)
        yield i0, i1, box_flat_indices(offsets[i0:i1], shapes[i0:i1], lo[i0:i1], size[i0:i1])
        i0 = i1


# --------------------------------------------------------------------------- spatial index


class BoxPairs(NamedTuple):
    """Result of ``BinIndex.pairs``: one row per intersecting (query box, indexed box, shift).

    a     : (m,) int64  row of the query table
    b     : (m,) int64  row of the indexed table
    shift : (m,) int64  row of the ``shifts`` array used (0 when no shifts were given)
    lo,hi : (m, dim) int64 the intersection box ``(query[a] + shifts[shift]) ∩ indexed[b]``,
            in the coordinates of the INDEXED table (subtract the shift to go back to the
            query's coordinates)
    Rows are sorted by ``(a, b, shift)``.
    """

    a: np.ndarray
    b: np.ndarray
    shift: np.ndarray
    lo: np.ndarray
    hi: np.ndarray


def _empty_pairs(dim):
    z = np.zeros(0, dtype=np.int64)
    zz = np.zeros((0, dim), dtype=np.int64)
    return BoxPairs(z, z.copy(), z.copy(), zz, zz.copy())


class BinIndex:
    """Uniform-bin spatial index over one box table (closed integer boxes).

    The integer space is cut in bins of ``width`` points per axis (default: the mean box
    extent, rounded up). Every box is registered in each bin it overlaps; the (bin, box)
    entries are kept sorted by a linear bin key, so only occupied bins cost memory.
    Build O(N) entries (N log N sort). A query visits the bins of the query box only;
    a pair found in several bins is reported once (in the bin holding the lower corner
    of the intersection). Query boxes covering more than ``big_query_bins`` bins are
    answered by a vectorised scan of the whole table instead.
    """

    big_query_bins = 4096
    block = 1 << 16  # query boxes processed per numpy pass, bounds temporary memory

    def __init__(self, lower, upper, width=None):
        self.lower, self.upper = as_box_table(lower, upper)
        self.n, self.dim = self.lower.shape
        if self.n == 0:
            self.width = np.ones(self.dim, dtype=np.int64)
            self.keys = np.zeros(0, dtype=np.int64)
            self.ids = np.zeros(0, dtype=np.int64)
            self.starts = None
            return
        if width is None:
            ext = (self.upper - self.lower + 1).mean(axis=0)
            width = np.maximum(1, np.ceil(ext)).astype(np.int64)
        self.width = np.broadcast_to(np.asarray(width, dtype=np.int64), (self.dim,)).copy()
        self.bbox_lo = self.lower.min(axis=0)
        self.bbox_hi = self.upper.max(axis=0)
        self.bmin = self.bbox_lo // self.width
        self.bmax = self.bbox_hi // self.width
        nbins = self.bmax - self.bmin + 1
        if math.prod(int(v) for v in nbins) >= 2**62:
            raise ValueError("BinIndex: too many bins for an int64 key; pass a larger width")
        self.mult = np.ones(self.dim, dtype=np.int64)
        if self.dim > 1:
            self.mult[1:] = np.cumprod(nbins[:-1])
        self._lower_t = self.lower.T.copy()  # column copies for the pair kernel
        self._upper_t = self.upper.T.copy()
        ids, keys = self._entries(self.lower, self.upper)
        order = np.argsort(keys, kind="stable")
        self.keys = keys[order]
        self.ids = ids[order]
        # Dense table of where each bin's entries start, when it is small compared with the
        # entries: a lookup is then two gathers instead of two binary searches (several times
        # faster on 10^5-10^6 queries). Sparse layouts (few boxes in a huge bounding box) keep
        # the binary search.
        total_bins = math.prod(int(v) for v in nbins)
        self.starts = None
        if total_bins <= max(4 * len(self.keys), 2**20):
            starts = np.zeros(total_bins + 1, dtype=np.int64)
            np.cumsum(np.bincount(self.keys, minlength=total_bins), out=starts[1:])
            self.starts = starts

    # ---- helpers
    def _bin_ranges(self, lo, hi):
        """Bin ranges of boxes, clipped to the occupied bin range; count 0 if outside."""
        blo = np.maximum(lo // self.width, self.bmin)
        bhi = np.minimum(hi // self.width, self.bmax)
        nb = np.clip(bhi - blo + 1, 0, None)
        return blo, nb

    def _entries(self, lo, hi, rows=None):
        """(row, bin key) of every bin overlapped by every box (row = position in lo)."""
        blo, nb = self._bin_ranges(lo, hi)
        counts = np.prod(nb, axis=1)
        row = np.repeat(np.arange(len(lo), dtype=np.int64), counts)
        r = concat_ranges(np.zeros(len(lo), dtype=np.int64), counts)
        key = np.zeros(len(row), dtype=np.int64)
        for a in range(self.dim):
            s = nb[row, a]
            key += (blo[row, a] - self.bmin[a] + r % s) * self.mult[a]
            r //= s
        return (row if rows is None else rows[row]), key

    def _key_of_point(self, pts):
        return ((pts // self.width - self.bmin) * self.mult).sum(axis=1)

    @property
    def nbytes(self):
        """Bytes held by the index tables."""
        n = self.keys.nbytes + self.ids.nbytes + 2 * (self.lower.nbytes + self.upper.nbytes)
        return int(n + (0 if self.starts is None else self.starts.nbytes))

    # ---- queries
    def query(self, lo, hi):
        """Sorted int64 rows of the indexed boxes that share a point with the box ``[lo, hi]``."""
        lo = np.asarray(lo, dtype=np.int64).reshape(1, self.dim)
        hi = np.asarray(hi, dtype=np.int64).reshape(1, self.dim)
        return np.unique(self.pairs(lo, hi).b)

    def pairs(self, lo, hi, shifts=None, sort=True):
        """All intersecting pairs between a query table and the indexed table.

        lo, hi : (m, dim) int64 query boxes (closed).
        shifts : optional (k, dim) int64 integer translations (e.g. ``periodic_shifts``);
                 every query box is tested at every shift: ``query + shift`` vs indexed.
                 Without shifts, a single zero shift is used.
        sort   : True (default): rows sorted by ``(a, b, shift)``. False: same rows in an
                 unspecified order (saves a sort of all rows when the caller does not care).
        Returns ``BoxPairs`` (see there). A box intersecting its own periodic image is
        reported (when the table indexed is the query table and the shift is non-zero);
        with a zero shift a box of the same table pairs with itself too — drop ``a == b``
        rows if that is not wanted.
        """
        qlo, qhi = as_box_table(lo, hi)
        if qlo.shape[1] != self.dim:
            raise ValueError(f"query boxes have dim {qlo.shape[1]}, index has dim {self.dim}")
        if shifts is None:
            shifts = np.zeros((1, self.dim), dtype=np.int64)
        shifts = np.atleast_2d(np.asarray(shifts, dtype=np.int64))
        if self.n == 0 or len(qlo) == 0:
            return _empty_pairs(self.dim)
        reach = self._reaching_rows(qlo, qhi, shifts)
        parts = []
        for k, s in enumerate(shifts):
            near = reach(s)  # only query boxes that reach the indexed bounding box can match
            slo, shi = qlo[near] + s, qhi[near] + s
            for b0 in range(0, len(near), self.block):
                rows = near[b0 : b0 + self.block]
                blk = slice(b0, b0 + self.block)
                a, b, lo_, hi_ = self._pairs_block(slo[blk], shi[blk], rows)
                parts.append((a, b, np.full(len(a), k, dtype=np.int64), lo_, hi_))
        if not parts:
            return _empty_pairs(self.dim)
        a, b, sh, ilo, ihi = (np.concatenate(x) for x in zip(*parts))
        if not sort:
            return BoxPairs(a, b, sh, ilo, ihi)
        order = np.lexsort((sh, b, a))
        return BoxPairs(a[order], b[order], sh[order], ilo[order], ihi[order])

    def _reaching_rows(self, qlo, qhi, shifts):
        """Function ``shift -> sorted rows`` of the query boxes that, moved by ``shift``, meet
        the indexed bounding box.

        Why: with periodic shifts, a shifted query box can only meet the indexed boxes when it
        sits near the domain boundary along every shifted axis. Testing each axis and shift
        value once, then intersecting the (small) sets, avoids a pass over every query box for
        each of the 3**dim shifts.
        """
        per_axis = {}

        def rows_for(a, v):
            if (a, v) not in per_axis:
                ok = (qhi[:, a] + v >= self.bbox_lo[a]) & (qlo[:, a] + v <= self.bbox_hi[a])
                per_axis[(a, v)] = np.nonzero(ok)[0]
            return per_axis[(a, v)]

        def reach(s):
            sets = [rows_for(a, int(s[a])) for a in range(self.dim)]
            i = min(range(self.dim), key=lambda a: len(sets[a]))
            rows = sets[i]
            for a in range(self.dim):
                if a != i and len(rows):
                    v = int(s[a])
                    keep = (qhi[rows, a] + v >= self.bbox_lo[a]) & (
                        qlo[rows, a] + v <= self.bbox_hi[a])
                    rows = rows[keep]
            return rows

        return reach

    def _pairs_block(self, qlo, qhi, rows):
        _, nb = self._bin_ranges(qlo, qhi)
        nbins = np.prod(nb, axis=1)
        big = nbins > self.big_query_bins
        out = []
        if (~big).any():
            small = np.nonzero(~big)[0]
            qrow, key = self._entries(qlo[small], qhi[small], rows=small)
            if self.starts is not None:
                s, cnt = self.starts[key], self.starts[key + 1] - self.starts[key]
            else:
                s = np.searchsorted(self.keys, key, side="left")
                cnt = np.searchsorted(self.keys, key, side="right") - s
            a = np.repeat(qrow, cnt)
            kb = np.repeat(key, cnt)
            b = self.ids[concat_ranges(s, cnt)]
            # Intersect axis by axis on 1D columns: gathers of contiguous columns are about
            # 4x faster than gathers of (m, dim) rows on the millions of candidates.
            qlo_t, qhi_t = qlo.T.copy(), qhi.T.copy()
            ok = np.ones(len(a), dtype=bool)
            corner_key = np.zeros(len(a), dtype=np.int64)
            los, his = [], []
            for d in range(self.dim):
                lo_d = np.maximum(qlo_t[d][a], self._lower_t[d][b])
                hi_d = np.minimum(qhi_t[d][a], self._upper_t[d][b])
                ok &= lo_d <= hi_d
                corner_key += (lo_d // self.width[d] - self.bmin[d]) * self.mult[d]
                los.append(lo_d)
                his.append(hi_d)
            # report a pair only from the bin that holds the intersection's lower corner
            ok &= corner_key == kb
            sel = np.nonzero(ok)[0]
            out.append((a[sel], b[sel], np.stack([v[sel] for v in los], axis=1),
                        np.stack([v[sel] for v in his], axis=1)))
        for i in np.nonzero(big)[0]:  # rare huge query boxes: scan the whole table
            lo, hi, ok = intersect(qlo[i], qhi[i], self.lower, self.upper)
            b = np.nonzero(ok)[0]
            out.append((np.full(len(b), i, dtype=np.int64), b, lo[b], hi[b]))
        a, b, lo, hi = (np.concatenate(x) for x in zip(*out))
        return rows[a], b, lo, hi


def intersect_tables(alo, ahi, blo, bhi, shifts=None):
    """All intersecting pairs between box table A (queries) and box table B (indexed).

    Convenience for ``BinIndex(blo, bhi).pairs(alo, ahi, shifts)``; see ``BoxPairs``.
    """
    return BinIndex(blo, bhi).pairs(alo, ahi, shifts)


# --------------------------------------------------------------------------- physical boxes


_NEW_SYNTAX = (
    "boxes are physical and named: box=dict(x=(x0, x1), y=(y0, y1)) in the units of the "
    "simulation (the same as the node coordinates). Axes left out mean the whole domain. "
    "To select by level-0 cell indices as before, use selection_box=Box(...) (deprecated)."
)


def _snap(q):
    """Snap values within 1e-9 (relative) of an integer to that integer.

    Why: 0.6 / 0.2 == 2.9999999999999996 must give cell 3, not 2.
    """
    r = np.round(q)
    close = np.abs(q - r) <= 1e-9 * np.maximum(1.0, np.abs(q))
    return np.where(close, r, q)


@dataclass(frozen=True)
class PhysBox:
    """An axis-aligned box in physical coordinates, ``lo[a] <= x_a <= hi[a]``.

    lo, hi : tuples of ``dim`` floats; ``-inf`` / ``+inf`` for an axis the user left free.

    Build it with ``PhysBox.parse(dict(x=(x0, x1), y=(y0, y1)), dim)``. Conversion to the
    cells of one level (``cells``), for a level of mesh size ``dl`` and ``N`` cells:
    ``[floor(x0 / dl), ceil(x1 / dl) - 1]`` per axis, clipped to ``[0, N - 1]``; quotients within
    1e-9 of an integer are snapped to it first. The result has at least one cell: a degenerate
    box ``x0 == x1`` selects the cell containing that plane, or both cells when the plane is a
    cell face. A patch is selected when its cell box intersects that cell box, so a patch that
    only touches the box edge is not selected.

    Example: dl = 0.5, N = 10, ``x=(1.2, 2.0)`` -> cells ``[2, 3]`` (cell 2 spans [1.0, 1.5],
    cell 3 spans [1.5, 2.0]); a patch with cells ``[4, 7]`` (starting at x = 2.0) is not
    selected. ``x=(1.0, 1.0)`` -> cells ``[1, 2]``; ``x=(1.2, 1.2)`` -> cells ``[2, 2]``.
    """

    lo: tuple
    hi: tuple

    def __post_init__(self):
        lo = tuple(float(v) for v in self.lo)
        hi = tuple(float(v) for v in self.hi)
        if len(lo) != len(hi) or not 1 <= len(lo) <= 3:
            raise ValueError(f"PhysBox needs 1 to 3 axes, got lo={lo} hi={hi}")
        for a, (l, h) in enumerate(zip(lo, hi)):
            if math.isnan(l) or math.isnan(h) or l > h:
                raise ValueError(f"PhysBox: axis {AXES[a]} needs x0 <= x1, got ({l}, {h})")
        object.__setattr__(self, "lo", lo)
        object.__setattr__(self, "hi", hi)

    @property
    def dim(self):
        return len(self.lo)

    @classmethod
    def parse(cls, box, dim):
        """Turn the user's ``box=`` argument into a PhysBox (or None for "whole domain").

        Accepts None, a PhysBox, or a dict with keys among ``x``, ``y``, ``z`` (the first
        ``dim`` of them) mapping to ``(x0, x1)`` pairs of numbers; ``None`` as a bound means
        unbounded on that side. Raises TypeError for anything else (in particular the former
        ``Box`` objects and bare tuples), with the new syntax in the message; ValueError for
        unknown axes or ``x0 > x1``.
        """
        if box is None:
            return None
        if isinstance(box, PhysBox):
            if box.dim != dim:
                raise ValueError(f"box has {box.dim} axes, the data has {dim}")
            return box
        if not isinstance(box, dict):
            raise TypeError(f"box={box!r} is not accepted: {_NEW_SYNTAX}")
        allowed = AXES[:dim]
        unknown = [k for k in box if k not in allowed]
        if unknown:
            raise ValueError(f"box: unknown axis {unknown}; the data has axes {allowed}")
        lo, hi = [-math.inf] * dim, [math.inf] * dim
        for a, name in enumerate(allowed):
            if name not in box:
                continue
            bounds = box[name]
            if isinstance(bounds, (str, bytes)) or not hasattr(bounds, "__len__"):
                raise TypeError(
                    f"box[{name!r}]={bounds!r}: expected a pair (x0, x1). {_NEW_SYNTAX}"
                )
            if len(bounds) != 2:
                raise ValueError(f"box[{name!r}] must be a pair (x0, x1), got {bounds!r}")
            if bounds[0] is not None:
                lo[a] = float(bounds[0])
            if bounds[1] is not None:
                hi[a] = float(bounds[1])
        return cls(tuple(lo), tuple(hi))

    @classmethod
    def from_cells(cls, lower, upper, dl):
        """PhysBox spanning the cells ``[lower, upper]`` of a level of mesh size ``dl``."""
        dl = np.asarray(dl, dtype=np.float64)
        lo = np.asarray(lower, dtype=np.float64) * dl
        hi = (np.asarray(upper, dtype=np.float64) + 1) * dl
        return cls(tuple(lo), tuple(hi))

    def to_dict(self):
        """``dict(x=(x0, x1), ...)`` (free axes as ``(-inf, inf)``)."""
        return {AXES[a]: (self.lo[a], self.hi[a]) for a in range(self.dim)}

    def intersect(self, other):
        """Intersection with another PhysBox (or None = whole domain); None if empty."""
        if other is None:
            return self
        if other.dim != self.dim:
            raise ValueError(f"cannot intersect boxes of {self.dim} and {other.dim} axes")
        lo = tuple(max(a, b) for a, b in zip(self.lo, other.lo))
        hi = tuple(min(a, b) for a, b in zip(self.hi, other.hi))
        if any(l > h for l, h in zip(lo, hi)):
            return None
        return PhysBox(lo, hi)

    def clipped(self, domain):
        """Bounds as two float64 ``(dim,)`` arrays, free sides replaced by ``0`` / ``domain``."""
        domain = np.asarray(domain, dtype=np.float64)
        lo = np.maximum(np.asarray(self.lo), 0.0)
        hi = np.minimum(np.asarray(self.hi), domain)
        return lo, hi

    # DESIGN-NOTE: beyond DESIGN's rule: quotients are snapped to integers within 1e-9, the box
    # is clipped to the domain first, and a box entirely outside the domain raises ValueError.
    def cells(self, dl, domain_cells):
        """Inclusive cell box ``(lo, hi)`` int64 ``(dim,)`` of one level (rule in the class doc).

        dl : (dim,) mesh size of the level; domain_cells : (dim,) number of cells of the level.
        Raises ValueError when the box lies outside the domain on some axis.
        """
        dl = np.asarray(dl, dtype=np.float64)
        n = np.asarray(domain_cells, dtype=np.int64)
        lo, hi = self.clipped(n * dl)
        if (lo > hi).any():
            a = int(np.nonzero(lo > hi)[0][0])
            raise ValueError(
                f"box {self.to_dict()} lies outside the domain [0, {n[a] * dl[a]}] along "
                f"{AXES[a]}; choose bounds inside the domain"
            )
        q0, q1 = _snap(lo / dl), _snap(hi / dl)
        clo = np.floor(q0).astype(np.int64)
        chi = np.ceil(q1).astype(np.int64) - 1
        # degenerate box on a cell face: take both cells around the plane
        face = (chi < clo) & (q0 == np.floor(q0))
        clo = np.where(face, clo - 1, clo)
        chi = np.maximum(chi, clo + np.where(face, 1, 0))
        return np.clip(clo, 0, n - 1), np.clip(chi, 0, n - 1)

    def selects(self, lower, upper, dl, domain_cells):
        """Bool ``(n,)``: which patches (cell boxes ``lower``, ``upper``) the box selects."""
        clo, chi = self.cells(dl, domain_cells)
        return ((np.asarray(upper) >= clo) & (np.asarray(lower) <= chi)).all(axis=1)

    def __repr__(self):
        parts = ", ".join(f"{AXES[a]}=({self.lo[a]:g}, {self.hi[a]:g})" for a in range(self.dim))
        return f"PhysBox({parts})"
