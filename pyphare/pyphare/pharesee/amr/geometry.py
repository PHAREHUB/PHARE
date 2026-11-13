"""Patch layout of AMR levels: where every patch lives in the flat per-level buffers.

Index conventions (DESIGN 4.1): level ``l`` has cells ``0 .. N_l - 1`` per axis,
``N_l = N_0 * 2**l``; a patch is an inclusive CELL box ``[lower, upper]``; its nodes are
``lower .. upper + 1`` (node shape ``upper - lower + 2``); node ``i`` sits at ``i * dl_l``.

Storage (DESIGN 4.2): one flat buffer per (time, level, component). Patches are stored in
STORE order = stable sort of the file-order patches by node shape ``(nx, ny, nz)``
(lexicographic, x first); inside a patch, values are x fastest. ``file_index`` and
``store_index`` translate between store slots and file-order patch numbers.

Nothing here reads files and nothing caches an array with one entry per node.
"""

from __future__ import annotations

import hashlib
import os
import threading
import weakref
from collections import OrderedDict
from dataclasses import dataclass
from typing import NamedTuple

import numpy as np

from .boxes import (
    BinIndex,
    PhysBox,
    as_box_table,
    box_flat_indices,
    coarsen_inner,
    concat_ranges,
    fortran_strides,
    periodic_shifts,
)


def _readonly(a):
    a = np.ascontiguousarray(a)
    a.flags.writeable = False
    return a


# DESIGN-NOTE: the key also covers domain_cells and periodic (DESIGN lists level, boxes, dl):
# two layouts that differ only there must not share buffers or face tables.
def geometry_key(level, lower, upper, dl, domain_cells, periodic):
    """Digest (bytes) identifying a level layout; ``lower``/``upper`` in FILE order.

    Two levels with equal keys have identical patches in identical order, mesh size,
    domain and periodicity, so their flat buffers are interchangeable.
    """
    h = hashlib.blake2b(digest_size=20)
    lower = np.ascontiguousarray(lower, dtype=np.int64)
    upper = np.ascontiguousarray(upper, dtype=np.int64)
    h.update(np.array([level, lower.shape[0], lower.shape[1]], dtype=np.int64).tobytes())
    h.update(lower.tobytes())
    h.update(upper.tobytes())
    h.update(np.ascontiguousarray(dl, dtype=np.float64).tobytes())
    h.update(np.ascontiguousarray(domain_cells, dtype=np.int64).tobytes())
    h.update(np.ascontiguousarray(periodic, dtype=bool).tobytes())
    return h.digest()


class EndLayers(NamedTuple):
    """First and last node layer of every patch along one axis (see LevelGeometry.end_layers).

    lo, hi : (m,) int64 flat node indices; ``lo[k]`` is on the first layer (local index 0)
             and ``hi[k]`` on the last layer (local index ``n_axis - 1``) of the same node line.
    stride : (m,) int64 flat distance between consecutive nodes along the axis in that patch
             (``hi[k] = lo[k] + (count[k] - 1) * stride[k]``).
    count  : (m,) int64 number of nodes of the patch along the axis (>= 2).
    Order: patches in store order; inside a patch, the layer's nodes x fastest over the
    remaining axes.
    """

    lo: np.ndarray
    hi: np.ndarray
    stride: np.ndarray
    count: np.ndarray


# DESIGN-NOTE: DESIGN 7.2 writes d[face] = (f[outer] - f[inner]) / (2 dl), which is right on hi
# faces only (on a lo face outer is at -dl); the table carries ``sign`` (and ``plus``/``minus``).
# ``one_sided_lo/hi`` are read as one bool per patch for this axis. A neighbour must hold both
# the face node and the outer node (the patches touch at the face node); a patch holding only
# the outer node across a one-cell gap does not count.
@dataclass(frozen=True)
class FaceTable:
    """Same-level neighbours of the face nodes of every patch along one axis.

    A *face node* of a patch along ``axis`` is a node on its first (``lo`` face) or last
    (``hi`` face) node layer along that axis. The table lists every face node for which some
    patch of the same level (another patch, or a periodic image of any patch including the
    patch itself) holds both the face node and the node one step OUTSIDE the patch along
    ``axis`` — i.e. the edge between them lies in that patch's node box. For these nodes a
    centred difference is possible with same-level data only.

    Arrays (all of length ``m``, one row per listed face node, sorted by ``face``):
      face  : int64 flat index of the face node in the level's buffer (its own patch)
      outer : int64 flat index, in the NEIGHBOUR's patch, of the node one step outside:
              position ``x_face - dl`` on a lo face, ``x_face + dl`` on a hi face (modulo the
              period through a periodic boundary)
      sign  : int8, -1 on a lo face, +1 on a hi face
      inner : int64 flat index of the node one step INSIDE the face's patch
              (``face + stride`` on a lo face, ``face - stride`` on a hi face, i.e.
              ``face - sign * stride``). Not stored: computed on every access (one pass
              over the rows); ``inner_rows(a, b)`` computes a block of rows only.
      patch : int64 store slot of the face's patch. Not stored: computed on every access.
    Per patch (store order, length npatch):
      one_sided_lo, one_sided_hi : bool, True when at least one node of that face is NOT
              listed (it keeps the one-sided formula: domain boundary or coarse-fine face).
      stride    : int64 flat distance between neighbours along ``axis`` in that patch
      row_start : int64 (npatch + 1,): rows of patch ``p`` are ``row_start[p]:row_start[p+1]``
              (rows sorted by ``face`` are grouped by patch, in store order)

    Formulas, with ``h = dl[axis]``:
      first derivative   d[face] = sign * (f[outer] - f[inner]) / (2 h)
                         == (f[plus] - f[minus]) / (2 h)  with the ``plus``/``minus`` properties
      second derivative  d2[face] = (f[outer] - 2 f[face] + f[inner]) / h**2
    Nodes shared by two patches appear once per patch (two flat indices) and receive the
    same value from both rows (same two data values, possibly from different patches).

    Worked example (1D, dl = 1, periodic, N = 8): patches A cells [0, 3] (nodes 0..4, store
    flat 0..4) and B cells [4, 7] (nodes 4..8, flat 5..9). Rows (face, inner, outer, sign):
    (0, 1, 8, -1) A's lo node 0, outside neighbour x=-1 is B's node 7 through the period;
    (4, 3, 6, +1) A's hi node 4, outside x=5 is B's node 5 (flat 6);
    (5, 6, 3, -1) B's lo node 4, outside x=3 is A's node 3; (9, 8, 1, +1) B's hi node 8 ->
    x=9 wraps to A's node 1. Without periodicity only the two middle rows remain and
    ``one_sided_lo = [True, False]``, ``one_sided_hi = [False, True]``.
    """

    axis: int
    face: np.ndarray
    outer: np.ndarray
    sign: np.ndarray
    one_sided_lo: np.ndarray
    one_sided_hi: np.ndarray
    stride: np.ndarray
    row_start: np.ndarray

    def __len__(self):
        return len(self.face)

    def _row_stride(self, a, b):
        """Stride of the patch of every row in [a, b)."""
        p0 = int(np.searchsorted(self.row_start, a, side="right")) - 1
        p1 = int(np.searchsorted(self.row_start, b, side="left"))
        counts = np.minimum(self.row_start[p0 + 1 : p1 + 1], b) - np.maximum(
            self.row_start[p0:p1], a)
        return np.repeat(self.stride[p0:p1], counts)

    def inner_rows(self, a, b):
        """``inner[a:b]`` computed for these rows only (int64, one pass over b - a rows)."""
        a, b = max(0, int(a)), min(len(self.face), int(b))
        if b <= a:
            return np.zeros(0, dtype=np.int64)
        return self.face[a:b] - self.sign[a:b] * self._row_stride(a, b)

    @property
    def inner(self):
        """int64 flat index of the node one step inside the face's patch (computed)."""
        return self.inner_rows(0, len(self.face))

    @property
    def patch(self):
        """int64 store slot of the patch of every ``face`` node (computed, not stored)."""
        npatch = len(self.stride)
        return np.repeat(np.arange(npatch, dtype=np.int64), np.diff(self.row_start))

    @property
    def plus(self):
        """Flat index of the node at ``+dl`` along the axis (outer on hi faces, inner on lo)."""
        return np.where(self.sign > 0, self.outer, self.inner)

    @property
    def minus(self):
        """Flat index of the node at ``-dl`` along the axis (inner on hi faces, outer on lo)."""
        return np.where(self.sign > 0, self.inner, self.outer)

    @property
    def nbytes(self):
        """Bytes held by the table: 17 per row (face, outer, sign) + 18 per patch."""
        return int(
            sum(
                getattr(self, k).nbytes
                for k in ("face", "outer", "sign", "one_sided_lo", "one_sided_hi", "stride",
                          "row_start")
            )
        )


@dataclass(frozen=True)
class CoverTable:
    """Parts of the patches of level ``level`` covered by level ``level + 1``.

    One row per *piece* = intersection of one patch's cells with the cells of one finer
    patch coarsened by 2 (only coarse cells fully covered by that fine patch count; for
    ratio-2 aligned fine boxes, as PHARE writes them, that is the whole coarsened box).
    Pieces are disjoint: same-level patches do not share cells.

      slot      : (npieces,) int64 store slot of the coarse patch
      lo, hi    : (npieces, dim) int64 inclusive cell sub-box, LOCAL to the patch
                  (0-based: global cell = ``lo + geom.lower[slot]``)
      fine_slot : (npieces,) int64 store slot, in level ``level + 1``, of the covering patch
    Rows sorted by ``slot`` (then by fine_slot). The nodes of piece ``k`` are the local
    nodes ``lo[k] .. hi[k] + 1``; its flat node indices come from ``node_indices``.

    Worked example (1D): coarse patch slot 0 with cells [0, 9], fine patch cells [6, 11]
    (coarse cells 3..5): one piece ``slot=0, lo=[3], hi=[5]`` -> covered nodes 3..6 of the
    patch (local), i.e. the interval [3 dl, 6 dl].
    """

    level: int
    slot: np.ndarray
    lo: np.ndarray
    hi: np.ndarray
    fine_slot: np.ndarray

    def __len__(self):
        return len(self.slot)

    def node_indices(self, geom, closed=True):
        """Flat node indices (int64) of every piece in ``geom`` (the coarse level).

        closed=True: nodes ``lo .. hi + 1`` (the piece's nodes including its boundary);
        closed=False: nodes strictly inside (``lo + 1 .. hi``). Pieces one after the other,
        x fastest; a node on the border of two adjacent pieces appears twice.
        """
        if closed:
            lo, size = self.lo, self.hi - self.lo + 2
        else:
            lo, size = self.lo + 1, self.hi - self.lo
        return box_flat_indices(geom.offset[self.slot], geom.shape[self.slot], lo, size)

    @property
    def nbytes(self):
        return int(self.slot.nbytes + self.lo.nbytes + self.hi.nbytes + self.fine_slot.nbytes)


class LevelGeometry:
    """Patch layout of one AMR level (immutable, shared by every quantity/time with equal boxes).

    Attributes (all int64 tables are read-only, STORE order unless stated):
      level, dim        int
      lower, upper      (npatch, dim) inclusive cell boxes
      shape             (npatch, dim) node shapes ``upper - lower + 2``
      size              (npatch,) nodes per patch; offset (npatch + 1,) start of each patch in
                        the flat buffer (``offset[-1] == nnodes``)
      nnodes            int, length of the level's flat buffers
      file_index        (npatch,) store slot -> file-order patch number
      store_index       (npatch,) file-order patch number -> store slot
      dl                (dim,) float64 mesh size OF THIS LEVEL; domain_cells (dim,) int64 cells
                        of this level; periodic (dim,) bool
      groups            structured array, one row per shape group (patches of equal shape,
                        contiguous in store order): fields ``shape`` (dim,), ``p0, p1`` (patch
                        slots [p0, p1)), ``n0, n1`` (flat nodes [n0, n1))
    Build with ``from_file_order``. Equality and hash follow ``key()``.
    """

    def __init__(self, level, lower, upper, dl, domain_cells, periodic):
        lower, upper = as_box_table(lower, upper)
        self.level = int(level)
        self.dim = lower.shape[1]
        self.dl = _readonly(np.broadcast_to(np.asarray(dl, dtype=np.float64), (self.dim,)))
        self.domain_cells = _readonly(
            np.broadcast_to(np.asarray(domain_cells, dtype=np.int64), (self.dim,))
        )
        self.periodic = _readonly(np.broadcast_to(np.asarray(periodic, dtype=bool), (self.dim,)))
        if len(lower) and ((lower < 0).any() or (upper >= self.domain_cells).any()):
            bad = int(np.nonzero(((lower < 0) | (upper >= self.domain_cells)).any(axis=1))[0][0])
            raise ValueError(
                f"level {level}: patch {bad} (file order) cells {lower[bad]}..{upper[bad]} "
                f"outside the domain 0..{self.domain_cells - 1}"
            )
        self._key = geometry_key(level, lower, upper, self.dl, self.domain_cells, self.periodic)
        shape = upper - lower + 2
        # lexsort uses the LAST key as primary: pass (nz, ny, nx) so that nx is primary
        order = np.lexsort(shape.T[::-1]) if len(shape) else np.zeros(0, dtype=np.int64)
        self.file_index = _readonly(order.astype(np.int64))
        store_index = np.empty_like(self.file_index)
        store_index[self.file_index] = np.arange(len(order), dtype=np.int64)
        self.store_index = _readonly(store_index)
        self.lower = _readonly(lower[order])
        self.upper = _readonly(upper[order])
        self.shape = _readonly(shape[order])
        self.size = _readonly(np.prod(self.shape, axis=1))
        offset = np.zeros(len(order) + 1, dtype=np.int64)
        np.cumsum(self.size, out=offset[1:])
        self.offset = _readonly(offset)
        self.nnodes = int(offset[-1])
        self.groups = self._make_groups()
        self._lock = threading.Lock()
        self._bins = None
        self._runs = {}
        self._faces = {}
        self._owner = None  # weakref to the GeometryCache that holds this geometry
        self._evicted = False  # True once a GeometryCache dropped it: face tables not kept

    @property
    def npatch(self):
        return len(self.lower)

    @classmethod
    def from_file_order(cls, level, lower, upper, dl, domain_cells, periodic):
        """Build the layout from patch boxes in FILE order.

        level        : level number (0 = coarsest)
        lower, upper : (npatch, dim) int inclusive cell boxes, file order
        dl           : (dim,) mesh size of THIS level (``dl0 / 2**level``)
        domain_cells : (dim,) number of cells of THIS level (``N0 * 2**level``)
        periodic     : (dim,) bool
        Raises ValueError for malformed boxes or boxes outside the domain. Use
        ``for_level`` to pass level-0 values instead.
        """
        # DESIGN-NOTE: DESIGN 4.2 does not say whether dl/domain_cells are level-0 or
        # this level's values; we take this level's (they are the attributes stored).
        return cls(level, lower, upper, dl, domain_cells, periodic)

    @classmethod
    def for_level(cls, level, lower, upper, dl0, domain_cells0, periodic):
        """``from_file_order`` with level-0 mesh size and cell counts (ratio 2 per level)."""
        r = 2**level
        dl = np.asarray(dl0, dtype=np.float64) / r
        n = np.asarray(domain_cells0, dtype=np.int64) * r
        return cls(level, lower, upper, dl, n, periodic)

    def file_order_boxes(self):
        """``(lower, upper)`` in FILE order (fresh arrays)."""
        return self.lower[self.store_index], self.upper[self.store_index]

    # ---------------------------------------------------------------- identity
    def key(self):
        """Digest of (level, lower, upper in file order, dl, domain_cells, periodic)."""
        return self._key

    def __hash__(self):
        return hash(self._key)

    def __eq__(self, other):
        return isinstance(other, LevelGeometry) and other._key == self._key

    def __reduce__(self):
        lower, upper = self.file_order_boxes()
        return (
            LevelGeometry,
            (self.level, lower, upper, np.array(self.dl), np.array(self.domain_cells),
             np.array(self.periodic)),
        )

    def __repr__(self):
        return (
            f"LevelGeometry(level={self.level}, dim={self.dim}, npatch={self.npatch}, "
            f"nnodes={self.nnodes}, groups={len(self.groups)})"
        )

    def nbytes(self):
        """Bytes held by this object, including cached tables (bin index, runs, faces)."""
        n = sum(
            getattr(self, k).nbytes
            for k in ("lower", "upper", "shape", "size", "offset", "file_index", "store_index",
                      "groups")
        )
        if self._bins is not None:
            n += self._bins.nbytes
        n += sum(sum(a.nbytes for a in r) for r in self._runs.values())
        n += sum(f.nbytes for f in self._faces.values())
        return int(n)

    # ---------------------------------------------------------------- groups, runs
    def _make_groups(self):
        dt = np.dtype(
            [("shape", np.int64, (self.dim,)), ("p0", np.int64), ("p1", np.int64),
             ("n0", np.int64), ("n1", np.int64)]
        )
        if self.npatch == 0:
            return _readonly(np.zeros(0, dtype=dt))
        starts = self._change_points(self.shape)
        p0 = starts
        p1 = np.append(starts[1:], self.npatch)
        g = np.zeros(len(p0), dtype=dt)
        g["shape"] = self.shape[p0]
        g["p0"], g["p1"] = p0, p1
        g["n0"], g["n1"] = self.offset[p0], self.offset[p1]
        return _readonly(g)

    @staticmethod
    def _change_points(keys):
        """Indices where a row of ``keys`` (n, k) differs from the previous one (0 included)."""
        if len(keys) == 0:
            return np.zeros(0, dtype=np.int64)
        diff = (keys[1:] != keys[:-1]).any(axis=1)
        return np.concatenate([[0], np.nonzero(diff)[0] + 1]).astype(np.int64)

    def runs(self, axis):
        """Stride runs along ``axis``: ``(stride, n0, n1)``, three int64 arrays of length R.

        Run ``r`` covers the flat nodes ``[n0[r], n1[r])``, made of whole consecutive patches
        whose node shapes agree on the axes BEFORE ``axis`` (axis 0: the whole level, stride 1;
        axis 1: equal nx, stride nx; axis 2: equal (nx, ny), stride nx*ny). Inside a run, for
        every node q that is NOT on the first or last layer along ``axis`` of its patch,
        ``q - stride[r]`` and ``q + stride[r]`` are its neighbours at ``-dl`` and ``+dl`` in the
        same patch. Nodes on those end layers get wrong neighbours from a flat shift and must be
        overwritten (``end_nodes``, then the ``FaceTable``). Centred-difference kernel::

            for s, a, b in zip(*geom.runs(axis)):
                d[a + s:b - s] = (f[a + 2*s:b] - f[a:b - 2*s]) / (2 * dl[axis])

        Example (2D): patches of node shapes (3, 2), (3, 4), (5, 2) in store order, offsets
        0, 6, 18 (28 nodes): ``runs(1) = ([3, 5], [0, 18], [18, 28])``;
        ``runs(0) = ([1], [0], [28])``.
        Empty level: three empty arrays.
        """
        if not 0 <= axis < self.dim:
            raise ValueError(f"axis {axis} out of range for a {self.dim}D level")
        with self._lock:
            if axis not in self._runs:
                if self.npatch == 0:
                    z = np.zeros(0, dtype=np.int64)
                    r = (z, z.copy(), z.copy())
                elif axis == 0:
                    r = (np.ones(1, np.int64), np.zeros(1, np.int64),
                         np.array([self.nnodes], np.int64))
                else:
                    p0 = self._change_points(self.shape[:, :axis])
                    p1 = np.append(p0[1:], self.npatch)
                    stride = np.prod(self.shape[p0, :axis], axis=1)
                    r = (stride, self.offset[p0], self.offset[p1])
                self._runs[axis] = tuple(_readonly(a) for a in r)
            return self._runs[axis]

    def end_layers(self, axis):
        """First/last node layer of every patch along ``axis`` as ``EndLayers``.

        Computed on request (one entry per node of those layers), not cached.
        One-sided 3-point kernel for the first layer (count >= 3)::

            lo, hi, s, n = geom.end_layers(axis)
            d[lo] = (-3 f[lo] + 4 f[lo + s] - f[lo + 2 s]) / (2 h)
            d[hi] = ( 3 f[hi] - 4 f[hi - s] + f[hi - 2 s]) / (2 h)
        and where ``n == 2``: ``d[lo] = d[hi] = (f[hi] - f[lo]) / h``.
        """
        if not 0 <= axis < self.dim:
            raise ValueError(f"axis {axis} out of range for a {self.dim}D level")
        lo_local = np.zeros_like(self.shape)
        size = self.shape.copy()
        size[:, axis] = 1
        lo = box_flat_indices(self.offset[:-1], self.shape, lo_local, size)
        per_patch = np.prod(size, axis=1)
        stride = np.repeat(fortran_strides(self.shape)[:, axis], per_patch)
        count = np.repeat(self.shape[:, axis], per_patch)
        return EndLayers(lo, lo + (count - 1) * stride, stride, count)

    # DESIGN-NOTE: end_nodes keeps DESIGN's (lo, hi) signature; kernels also need the stride
    # and node count per entry, given by end_layers.
    def end_nodes(self, axis):
        """``(lo, hi)`` int64 flat indices of the first/last node layer along ``axis``.

        Same as the first two fields of ``end_layers(axis)`` (see there for order and kernel).
        """
        e = self.end_layers(axis)
        return e.lo, e.hi

    # ---------------------------------------------------------------- views
    def view(self, flat, slot):
        """Patch ``slot`` (store order) of a flat buffer, indexed ``[x, y, z]``, no copy."""
        slot = int(slot)
        return flat[self.offset[slot] : self.offset[slot + 1]].reshape(
            tuple(self.shape[slot]), order="F"
        )

    def stack(self, flat, g):
        """Group ``g`` of a flat buffer as one array indexed ``[x, y, z, p]``, no copy.

        ``p`` runs over the group's patches in store order: ``stack(flat, g)[..., k]`` is
        ``view(flat, groups["p0"][g] + k)``.
        """
        row = self.groups[int(g)]
        n = int(row["p1"] - row["p0"])
        return flat[row["n0"] : row["n1"]].reshape(tuple(row["shape"]) + (n,), order="F")

    def coords(self, slot, axis):
        """Node coordinates (float64) of patch ``slot`` along ``axis``: ``(lower + i) * dl``."""
        slot = int(slot)
        i = np.arange(self.shape[slot, axis], dtype=np.float64)
        return (self.lower[slot, axis] + i) * self.dl[axis]

    def node_coords(self, axis):
        """Coordinate along ``axis`` of every node of the flat buffer, float64 ``(nnodes,)``.

        Computed on request (never cached): ~3 temporary int64 arrays of nnodes.
        """
        if self.npatch == 0:
            return np.zeros(0)
        patch = np.repeat(np.arange(self.npatch), self.size)
        local = np.arange(self.nnodes, dtype=np.int64) - self.offset[:-1][patch]
        stride = fortran_strides(self.shape)[:, axis]
        i = (local // stride[patch]) % self.shape[patch, axis]
        return (self.lower[patch, axis] + i) * self.dl[axis]

    # ---------------------------------------------------------------- spatial queries
    def bin_index(self):
        """BinIndex over the cell boxes (store order), built on first use and kept."""
        with self._lock:
            if self._bins is None:
                self._bins = BinIndex(self.lower, self.upper)
            return self._bins

    def patches_touching(self, lo, hi):
        """Sorted int64 store slots of patches whose cell box shares a cell with ``[lo, hi]``.

        lo, hi : (dim,) inclusive CELL indices of this level.
        """
        return self.bin_index().query(lo, hi)

    def patches_in(self, box):
        """Store slots of the patches selected by a ``PhysBox`` (None = all patches)."""
        if box is None:
            return np.arange(self.npatch, dtype=np.int64)
        lo, hi = box.cells(self.dl, self.domain_cells)
        return self.patches_touching(lo, hi)

    def faces(self, axis):
        """``FaceTable`` along ``axis`` (built on request, then cached under the cache cap).

        Cached only when this geometry plus ALL its face tables (one per axis, the missing
        ones estimated at this table's size) fit under the cap of the owning GeometryCache
        (else the default cache's cap). A geometry that a GeometryCache evicted keeps no
        table at all (nothing would account for it): each call builds the table again.
        Otherwise the table is returned without being kept
        and nothing else is evicted: derivatives call ``faces`` for every axis in turn, and
        keeping one table at the expense of the others would rebuild them on every call while
        flushing every other cached geometry. Sizes: ~17 bytes per listed face node, e.g.
        ~59 MiB per axis for 2*10^5 2D patches of 8x8 cells, ~165 MiB per axis for 10^5 3D
        patches of 6^3 cells.
        """
        if not 0 <= axis < self.dim:
            raise ValueError(f"axis {axis} out of range for a {self.dim}D level")
        with self._lock:
            table = self._faces.get(axis)
        if table is not None:
            return table
        table = build_faces(self, axis)
        owner = self._owner() if self._owner is not None else None
        if owner is None and self._evicted:
            return table
        cap = owner.max_bytes if owner is not None else default_cache().max_bytes
        with self._lock:
            cached = {a: t for a, t in self._faces.items() if a != axis}
            missing = self.dim - len(cached)  # this axis included
            whole = self.nbytes() + table.nbytes * missing
            keep = whole <= cap
            if keep:
                self._faces[axis] = table
        if keep and owner is not None:
            owner.trim(keep=self)
        return table

    # DESIGN-NOTE: coverage tables (CoverTable) are cached on HierGeometry objects, which
    # belong to fields, not on LevelGeometry: they are O(patches) and go with the field.
    def drop_cached_tables(self):
        """Forget the cached face tables (called when a cache evicts this geometry)."""
        with self._lock:
            self._faces.clear()


# -------------------------------------------------------------------- face tables


def build_faces(geom, axis):
    """Build the ``FaceTable`` of ``geom`` along ``axis`` (vectorised, O(face nodes)).

    Method: patch q can serve as neighbour of a node g along ``axis`` when q's node box holds
    g and g + e_axis, i.e. g lies in q's "edge box" (cells along ``axis``, nodes along the
    other axes). Hi faces query that table with their face nodes, lo faces with the nodes one
    step outside; periodic images through ``periodic_shifts`` (only boxes near the domain
    boundary are tried at a non-zero shift). Pairs come from ``BinIndex``.
    A face node with several neighbours (edges, corners) keeps the neighbour of smallest
    (store slot, shift row): the table does not depend on the order pairs are found in.
    """
    n, dim = geom.npatch, geom.dim
    empty = np.zeros(0, dtype=np.int64)
    if n == 0:
        f = np.zeros(0, dtype=bool)
        return FaceTable(axis, empty, empty.copy(), np.zeros(0, np.int8), f, f, empty.copy(),
                         np.zeros(1, np.int64))
    lower, upper, shape = geom.lower, geom.upper, geom.shape
    offset = geom.offset[:-1]
    strides = fortran_strides(shape)

    edge_lo, edge_hi = lower.copy(), upper + 1
    edge_hi[:, axis] = upper[:, axis]
    index = BinIndex(edge_lo, edge_hi)

    # query rows 0..n-1: hi faces (face nodes); rows n..2n-1: lo faces (outer nodes)
    q_lo = np.concatenate([lower, lower])
    q_hi = np.concatenate([upper + 1, upper + 1])
    q_lo[:n, axis] = q_hi[:n, axis] = upper[:, axis] + 1
    q_lo[n:, axis] = q_hi[n:, axis] = lower[:, axis] - 1
    shifts = periodic_shifts(geom.domain_cells, geom.periodic)
    pairs = index.pairs(q_lo, q_hi, shifts, sort=False)

    hi_side = pairs.a < n
    p = np.where(hi_side, pairs.a, pairs.a - n)
    q = pairs.b
    step = hi_side.astype(np.int64)
    # matched nodes in q coordinates are pairs.lo..pairs.hi: face nodes (hi side) or
    # outer nodes (lo side); local first node of the sub-box in q (outer) and in p (face)
    outer_lo = pairs.lo - lower[q]
    outer_lo[:, axis] += step
    face_lo = pairs.lo - shifts[pairs.shift] - lower[p]
    face_lo[:, axis] += 1 - step
    size = pairs.hi - pairs.lo + 1
    face, outer, row = _twin_flat_indices(
        offset[p] + (face_lo * strides[p]).sum(axis=1), strides[p],
        offset[q] + (outer_lo * strides[q]).sum(axis=1), strides[q], size,
    )
    pos = _first_by_face(face, q * len(shifts) + pairs.shift, row, geom.nnodes)
    face, outer, row = face[pos], outer[pos], row[pos]
    hi_node = hi_side[row]
    patch = p[row]
    sign = np.where(hi_node, 1, -1).astype(np.int8)
    # rows are sorted by face, hence grouped by patch in store order: ``inner`` is not
    # stored but recomputed from face, sign and the patch stride (FaceTable.inner_rows)
    row_start = np.zeros(n + 1, dtype=np.int64)
    np.cumsum(np.bincount(patch, minlength=n), out=row_start[1:])

    nodes_per_face = np.prod(np.delete(shape, axis, axis=1), axis=1) if dim > 1 else np.ones(n)
    listed_hi = np.bincount(patch[hi_node], minlength=n)
    listed_lo = np.bincount(patch[~hi_node], minlength=n)
    return FaceTable(
        axis,
        _readonly(face),
        _readonly(outer),
        _readonly(sign),
        _readonly(listed_lo < nodes_per_face),
        _readonly(listed_hi < nodes_per_face),
        _readonly(strides[:, axis]),
        _readonly(row_start),
    )


def _twin_flat_indices(base1, strides1, base2, strides2, size):
    """Flat indices of the same local sub-box enumeration in two patches at once.

    Row ``i``: a sub-box of ``size[i]`` nodes whose first node is at flat index ``base1[i]``
    in a patch of strides ``strides1[i]`` and at ``base2[i]`` in a patch of strides
    ``strides2[i]``. Returns ``(idx1, idx2, row)``: matching flat indices (node k of the
    sub-box in both patches) and the row of every node. Enumeration order: the axis with the
    most nodes in total fastest (one ``box_flat_indices`` pass shared by both patches).
    """
    n, dim = size.shape
    run = int(np.argmax(size.sum(axis=0)))
    others = [a for a in range(dim) if a != run]
    nlines = np.prod(size[:, others], axis=1) if others else np.ones(n, dtype=np.int64)
    nlines = np.where(size[:, run] > 0, nlines, 0)
    piece = np.repeat(np.arange(n, dtype=np.int64), nlines)
    k = concat_ranges(np.zeros(n, dtype=np.int64), nlines)
    s1, s2 = base1[piece], base2[piece]
    for a in others:
        m = size[piece, a]
        j = k % m
        s1 += j * strides1[piece, a]
        s2 += j * strides2[piece, a]
        k //= m
    length = size[piece, run]
    j = concat_ranges(np.zeros(len(piece), dtype=np.int64), length)
    idx1 = np.repeat(s1, length) + j * np.repeat(strides1[piece, run], length)
    idx2 = np.repeat(s2, length) + j * np.repeat(strides2[piece, run], length)
    return idx1, idx2, np.repeat(piece, length)


def _first_by_face(face, rank, row, nnodes):
    """Positions that sort ``face`` and keep one entry per distinct face node.

    Among entries of the same face node, the one whose ``rank[row]`` is smallest is kept.
    Why the packed sort: ``np.sort`` of int64 is several times faster than ``argsort``; the
    position rides in the low bits of ``face << bits | position`` when that fits in int64.
    """
    m = len(face)
    if m == 0:
        return np.zeros(0, dtype=np.int64)
    bits = int(m - 1).bit_length()
    if int(nnodes) << bits < 2**62:
        packed = np.sort((face << bits) | np.arange(m, dtype=np.int64))
        sorted_face = packed >> bits
        pos = packed & ((1 << bits) - 1)
    else:
        pos = np.argsort(face, kind="stable")
        sorted_face = face[pos]
    first = np.ones(m, dtype=bool)
    first[1:] = sorted_face[1:] != sorted_face[:-1]
    if not first.all():
        # face nodes with several neighbours: keep the entry of smallest rank
        multi = ~first
        multi[:-1] |= ~first[1:]
        # (ranks are distinct inside a group: one row per (neighbour, shift) and face node)
        idx = np.nonzero(multi)[0]  # entries of the multi-entry groups, group after group
        new_group = np.ones(len(idx), dtype=bool)
        new_group[1:] = first[idx[1:]]
        starts = np.nonzero(new_group)[0]
        r = rank[row[pos[idx]]]
        best = np.minimum.reduceat(r, starts)
        winner = r == np.repeat(best, np.diff(np.append(starts, len(idx))))
        pos[idx[starts]] = pos[idx[winner]]
    return pos[first]


# -------------------------------------------------------------------- hierarchy


class HierGeometry:
    """The levels present at ONE time, plus the extent of the field.

    levels : tuple of LevelGeometry, increasing level number; levels without patches are
             absent; level 0 must be present.
    extent : PhysBox, or None for the whole domain.
    """

    def __init__(self, levels, extent=None):
        levels = tuple(sorted(levels, key=lambda g: g.level))
        if not levels or levels[0].level != 0:
            raise ValueError("HierGeometry needs level 0")
        if len({g.level for g in levels}) != len(levels):
            raise ValueError("HierGeometry: a level number appears twice")
        dims = {g.dim for g in levels}
        if len(dims) != 1:
            raise ValueError(f"HierGeometry: levels of different dimensions {dims}")
        self.levels = levels
        self.dim = levels[0].dim
        if extent is not None and not isinstance(extent, PhysBox):
            raise TypeError("extent must be a PhysBox or None (use PhysBox.parse)")
        self.extent = extent
        self._covered = {}
        self._lock = threading.Lock()

    @property
    def dl0(self):
        return self.levels[0].dl

    @property
    def domain_cells0(self):
        return self.levels[0].domain_cells

    @property
    def periodic(self):
        return self.levels[0].periodic

    @property
    def domain(self):
        """Physical size of the domain, float64 ``(dim,)``."""
        return self.levels[0].domain_cells * self.levels[0].dl

    @property
    def finest_level(self):
        return self.levels[-1].level

    def level_numbers(self):
        return tuple(g.level for g in self.levels)

    def level(self, ilvl):
        """LevelGeometry of level number ``ilvl``; KeyError if that level has no patch."""
        for g in self.levels:
            if g.level == ilvl:
                return g
        raise KeyError(f"no level {ilvl}; levels present: {self.level_numbers()}")

    def has_level(self, ilvl):
        return any(g.level == ilvl for g in self.levels)

    def with_extent(self, extent):
        """Same levels, another extent (shares LevelGeometry objects and coverage tables)."""
        h = HierGeometry(self.levels, extent)
        h._covered = self._covered
        return h

    def key(self):
        """Tuple of level keys and the extent."""
        return tuple(g.key() for g in self.levels) + (self.extent,)

    def covered(self, ilvl):
        """``CoverTable`` of level ``ilvl``: its cells covered by level ``ilvl + 1``.

        Cost O(patches of both levels) through ``BinIndex``; cached on this object. Empty
        table for the finest level. Proper nesting is not required: only coarse cells fully
        covered by a fine patch are listed, so pieces stay disjoint for any non-overlapping
        fine patches.
        """
        with self._lock:
            if ilvl in self._covered:
                return self._covered[ilvl]
        coarse = self.level(ilvl)
        # DESIGN-NOTE: coarse cells only partly covered by one fine patch (fine boxes not
        # ratio-2 aligned, never written by PHARE) are not listed, which keeps pieces disjoint.
        z = np.zeros(0, dtype=np.int64)
        zz = np.zeros((0, self.dim), dtype=np.int64)
        if not self.has_level(ilvl + 1):
            table = CoverTable(ilvl, z, zz, zz.copy(), z.copy())
        else:
            fine = self.level(ilvl + 1)
            clo, chi = coarsen_inner(fine.lower, fine.upper)
            keep = np.nonzero((clo <= chi).all(axis=1))[0]
            pairs = coarse.bin_index().pairs(clo[keep], chi[keep])
            slot = pairs.b
            order = np.lexsort((keep[pairs.a], slot))
            slot = slot[order]
            table = CoverTable(
                ilvl,
                _readonly(slot),
                _readonly(pairs.lo[order] - coarse.lower[slot]),
                _readonly(pairs.hi[order] - coarse.lower[slot]),
                _readonly(keep[pairs.a][order]),
            )
        with self._lock:
            self._covered[ilvl] = table
        return table

    def __repr__(self):
        lv = ", ".join(f"L{g.level}:{g.npatch}p" for g in self.levels)
        return f"HierGeometry(dim={self.dim}, levels=[{lv}], extent={self.extent})"


# -------------------------------------------------------------------- cache


# DESIGN-NOTE: DESIGN 4.2 says 512 MiB. Measured face tables (17 B per row since ``inner`` is
# derived): ~59 MiB per axis at 2*10^5 2D patches, ~165 MiB per axis at 10^5 3D patches
# (~1.6 GiB per axis at 10^6 3D patches). 1 GiB
# keeps the three tables of a 10^5-patch 3D level; at 10^6 3D patches no sensible default
# can, and ``LevelGeometry.faces`` then returns them uncached instead of thrashing.
DEFAULT_CACHE_MB = 1024


def _env_cache_bytes():
    value = os.environ.get("PHARESEE_GEOM_CACHE_MB")
    if value is None or value == "":
        return DEFAULT_CACHE_MB * 2**20
    try:
        return int(float(value) * 2**20)
    except ValueError:
        raise ValueError(
            f"PHARESEE_GEOM_CACHE_MB={value!r} is not a number of MiB; unset it or give e.g. 512"
        ) from None


class GeometryCache:
    """Shared LevelGeometry objects, keyed by ``geometry_key``, LRU, capped in bytes.

    max_bytes : cap on the sum of ``LevelGeometry.nbytes()`` (default: env
                ``PHARESEE_GEOM_CACHE_MB`` in MiB, else ``DEFAULT_CACHE_MB`` = 1024 MiB). Face
                tables built on a cached geometry count toward the cap; evicting a geometry
                drops its face tables (objects still used elsewhere stay valid and rebuild them
                on request). A level whose face tables cannot all fit gets them uncached.
    Thread-safe.
    """

    def __init__(self, max_bytes=None):
        self.max_bytes = _env_cache_bytes() if max_bytes is None else int(max_bytes)
        self._items = OrderedDict()
        self._lock = threading.RLock()

    def get(self, level, lower, upper, dl, domain_cells, periodic):
        """The LevelGeometry for these FILE-order boxes (this level's dl/domain), built once."""
        key = geometry_key(level, lower, upper, dl, domain_cells, periodic)
        with self._lock:
            geom = self._items.get(key)
            if geom is not None:
                self._items.move_to_end(key)
                return geom
        geom = LevelGeometry.from_file_order(level, lower, upper, dl, domain_cells, periodic)
        return self.add(geom)

    def add(self, geom):
        """Insert ``geom`` (or return the equal one already cached)."""
        with self._lock:
            old = self._items.get(geom.key())
            if old is not None:
                self._items.move_to_end(geom.key())
                return old
            self._items[geom.key()] = geom
            if geom._owner is None or geom._owner() is None:
                geom._owner = weakref.ref(self)
                geom._evicted = False
            self.trim(keep=geom)
            return geom

    def nbytes(self):
        with self._lock:
            return sum(g.nbytes() for g in self._items.values())

    def trim(self, keep=None):
        """Evict least recently used geometries until the total fits ``max_bytes``."""
        with self._lock:
            total = self.nbytes()
            for key in list(self._items):
                if total <= self.max_bytes:
                    break
                geom = self._items[key]
                if geom is keep:
                    continue
                total -= geom.nbytes()
                del self._items[key]
                geom.drop_cached_tables()
                if geom._owner is not None and geom._owner() is self:
                    geom._owner = None
                    geom._evicted = True

    def clear(self):
        with self._lock:
            for geom in self._items.values():
                geom.drop_cached_tables()
                if geom._owner is not None and geom._owner() is self:
                    geom._owner = None
                    geom._evicted = True
            self._items.clear()

    def __len__(self):
        return len(self._items)

    def __contains__(self, key):
        return key in self._items


_default = {}
_default_lock = threading.Lock()


def default_cache():
    """The module-level GeometryCache (created on first use, reads PHARESEE_GEOM_CACHE_MB then)."""
    with _default_lock:
        if "cache" not in _default:
            _default["cache"] = GeometryCache()
        return _default["cache"]
