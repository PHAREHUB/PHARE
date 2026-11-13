"""Synthetic AMR layouts, analytic data on them, and an in-memory fake file backend.

For tests and benchmarks only. Everything is generated with numpy, without per-patch Python
loops, so that layouts of 10^5-10^6 patches are cheap to build.

* ``make_amr``      random but valid hierarchies (level 0 tiles the domain; finer levels are
                    ratio-2 aligned, non-overlapping, nested in the coarser level)
* ``make_large_amr`` many small patches for scaling benches
* ``fill``          in-memory LevelData from an analytic function of the node coordinates
* ``FakeBackend``   the Backend protocol on top of in-memory LevelData, counting the reads
"""

from __future__ import annotations

import itertools
import threading
import time as _time
import weakref

import numpy as np

from .boxes import PhysBox, concat_ranges, refine
from .errors import PhareseeIOError, TimeNotFoundError
from .geometry import HierGeometry, LevelGeometry
from .io import FileMeta, LevelIndex
from .store import LevelData

_DEFAULT_CELLS = {1: (64,), 2: (48, 40), 3: (16, 14, 12)}
_DEFAULT_WIDTH = {1: 8, 2: 6, 3: 4}


# --------------------------------------------------------------------------- tiling


def _split_axis(lo, hi, axis, width, rng, unit, jitter):
    """Cut every box along ``axis`` into segments of about ``width`` cells (vectorised).

    unit=2 keeps every cut on an even cell (ratio-2 alignment) given even box ends.
    jitter in [0, 0.25]: relative random displacement of the interior cuts.
    width may be a scalar or one value per box.
    """
    n = len(lo)
    ext = (hi[:, axis] - lo[:, axis] + 1) // unit  # extent in units
    w = np.maximum(np.broadcast_to(np.asarray(width) // unit, (n,)), 1)
    m = np.maximum(ext // w, 1)  # segments per box
    k = concat_ranges(np.zeros(n, np.int64), m + 1)  # 0..m per box
    box = np.repeat(np.arange(n), m + 1)
    cut = k * ext[box] // m[box]
    amp = (jitter * ext[box] / m[box]).astype(np.int64)
    interior = (k > 0) & (k < m[box])
    noise = rng.integers(-amp, amp + 1) if len(amp) else amp
    cut = np.where(interior, cut + noise, cut)
    # segments: consecutive cuts of the same box
    first = np.cumsum(m + 1) - (m + 1)
    seg_box = np.repeat(np.arange(n), m)
    seg_k = concat_ranges(first, m)  # index of the segment's left cut
    new_lo = lo[seg_box].copy()
    new_hi = hi[seg_box].copy()
    new_lo[:, axis] = lo[seg_box, axis] + cut[seg_k] * unit
    new_hi[:, axis] = lo[seg_box, axis] + cut[seg_k + 1] * unit - 1
    return new_lo, new_hi


def tile(lo, hi, width, rng, unit=1, jitter=0.25, singletons=False):
    """Tile every region box ``[lo, hi]`` (n, dim) into patches (x cut first, then y, z).

    Each column is cut independently along the next axis, so patches meet with offsets
    (T-junctions). Returns ``(lower, upper)`` int64. ``singletons=True`` draws a random target
    width per box at every axis so that most patch shapes are unique.
    """
    lo = np.atleast_2d(np.asarray(lo, np.int64))
    hi = np.atleast_2d(np.asarray(hi, np.int64))
    for axis in range(lo.shape[1]):
        w = width
        if singletons:
            w = rng.integers(max(unit, width // 2), width * 2 + 1, size=len(lo))
        lo, hi = _split_axis(lo, hi, axis, w, rng, unit, jitter)
    return lo, hi


def _peel_one_cell(lower, upper, rng, fraction=0.3):
    """Split a cell-wide slab off some patches so that 1-cell (2-node) patches exist."""
    dim = lower.shape[1]
    axis = rng.integers(0, dim, size=len(lower))
    rows = np.arange(len(lower))
    ext = upper[rows, axis] - lower[rows, axis] + 1
    pick = (rng.random(len(lower)) < fraction) & (ext >= 3)
    pick[np.argmax(ext)] |= ext.max() >= 3  # at least one
    idx = np.nonzero(pick)[0]
    a = axis[idx]
    slab_lo, slab_hi = lower[idx].copy(), upper[idx].copy()
    slab_lo[np.arange(len(idx)), a] = upper[idx, a]  # the last cell of the patch
    upper = upper.copy()
    upper[idx, a] -= 1
    return np.concatenate([lower, slab_lo]), np.concatenate([upper, slab_hi])


def _sub_regions(lo, hi, rng, fraction):
    """One or two random sub-boxes inside each region (level index space), disjoint."""
    out_lo, out_hi = [], []
    for rlo, rhi in zip(lo, hi):  # a handful of regions per level
        ext = rhi - rlo + 1
        parts = [(rlo, rhi)]
        ax = int(np.argmax(ext))
        if ext[ax] >= 4 and rng.random() < 0.5:
            mid = rlo[ax] + ext[ax] // 2
            a_hi, b_lo = rhi.copy(), rlo.copy()
            a_hi[ax], b_lo[ax] = mid - 1, mid
            parts = [(rlo, a_hi), (b_lo, rhi)]
        for plo, phi in parts:
            pext = phi - plo + 1
            size = np.maximum(1, np.round(pext * rng.uniform(fraction[0], fraction[1],
                                                             size=len(pext)))).astype(np.int64)
            start = plo + rng.integers(0, pext - size + 1)
            out_lo.append(start)
            out_hi.append(start + size - 1)
    return np.array(out_lo, np.int64), np.array(out_hi, np.int64)


def make_amr(dim, nlevels=2, seed=0, *, domain_cells=None, dl0=None, periodic=False,
             patch_cells=None, two_node=False, singletons=False, fraction=(0.3, 0.7),
             extent=None, shuffle=True):
    """Random valid AMR layout as a ``HierGeometry``.

    dim          1, 2 or 3
    nlevels      number of levels (level 0 included)
    seed         random seed (same seed -> same layout)
    domain_cells level-0 cells per axis (default 64 / 48x40 / 16x14x12)
    dl0          level-0 mesh size per axis (default 0.2 on every axis)
    periodic     bool or (dim,) bools
    patch_cells  target patch width in cells (default 8 / 6 / 4)
    two_node     True: some level-0 patches are 1 cell (2 nodes) wide along some axis
    singletons   True: most patch shapes are unique (groups of one patch)
    fraction     range of the relative size of each refined region inside its parent region
    extent       None or a box for ``PhysBox.parse`` (``dict(x=(x0, x1), ...)``)
    shuffle      True: file order is a random permutation (store order then differs)
    Level 0 tiles the whole domain; level l+1 tiles one or two random sub-regions of every
    level-l region, refined by 2, with even cuts (ratio-2 aligned, non-overlapping, nested).
    """
    rng = np.random.default_rng(seed)
    cells0 = np.asarray(domain_cells if domain_cells is not None else _DEFAULT_CELLS[dim],
                        np.int64)
    dl = np.broadcast_to(np.asarray(0.2 if dl0 is None else dl0, np.float64), (dim,))
    per = np.broadcast_to(np.asarray(periodic, bool), (dim,))
    width = patch_cells or _DEFAULT_WIDTH[dim]
    region_lo = np.zeros((1, dim), np.int64)
    region_hi = (cells0 - 1)[None, :]
    levels = []
    for ilvl in range(nlevels):
        unit = 1 if ilvl == 0 else 2
        lower, upper = tile(region_lo, region_hi, max(width, unit), rng, unit=unit,
                            singletons=singletons)
        if ilvl == 0 and two_node:
            lower, upper = _peel_one_cell(lower, upper, rng)
        if shuffle:
            perm = rng.permutation(len(lower))
            lower, upper = lower[perm], upper[perm]
        levels.append(LevelGeometry.for_level(ilvl, lower, upper, dl, cells0, per))
        if ilvl + 1 < nlevels:
            sub_lo, sub_hi = _sub_regions(region_lo, region_hi, rng, fraction)
            region_lo, region_hi = refine(sub_lo, sub_hi)
    return HierGeometry(levels, PhysBox.parse(extent, dim))


def make_large_amr(dim, npatch, nlevels=1, seed=0, width=8, periodic=False, dl0=0.2):
    """Hierarchy whose level 0 has about ``npatch`` patches of about ``width`` cells.

    Finer levels each refine one centred region of half the parent's extent per axis
    (their patch count is of the same order as level 0's in 2D). Vectorised: 2*10^5 patches
    in well under a second.
    """
    rng = np.random.default_rng(seed)
    per_axis = int(np.ceil(npatch ** (1.0 / dim)))
    cells0 = np.full(dim, per_axis * width, np.int64)
    dl = np.full(dim, dl0)
    region_lo = np.zeros((1, dim), np.int64)
    region_hi = (cells0 - 1)[None, :]
    levels = []
    for ilvl in range(nlevels):
        unit = 1 if ilvl == 0 else 2
        lower, upper = tile(region_lo, region_hi, width, rng, unit=unit)
        perm = rng.permutation(len(lower))
        levels.append(LevelGeometry.for_level(ilvl, lower[perm], upper[perm], dl, cells0,
                                              np.broadcast_to(periodic, (dim,))))
        ext = region_hi - region_lo + 1
        sub_lo = region_lo + ext // 4
        sub_hi = sub_lo + np.maximum(ext // 2, 1) - 1
        region_lo, region_hi = refine(sub_lo, sub_hi)
    return HierGeometry(levels)


# --------------------------------------------------------------------------- data


def fill(hier, f, names=None, dtype=np.float64, time=0.0):
    """In-memory LevelData (one per level) holding ``f`` at every node.

    f      callable ``f(x)`` / ``f(x, y)`` / ``f(x, y, z)`` on float64 coordinate arrays, or a
           dict ``{name: callable}`` for several components
    names  component names when ``f`` is a single callable (default ``("value",)``)
    Returns a list of LevelData in the order of ``hier.levels``.
    """
    funcs = dict(f) if isinstance(f, dict) else {(names or ("value",))[0]: f}
    out = []
    for geom in hier.levels:
        coords = [geom.node_coords(a) for a in range(geom.dim)]
        arrays = {
            name: np.asarray(fn(*coords), dtype=dtype).reshape(geom.nnodes).copy()
            for name, fn in funcs.items()
        }
        out.append(LevelData.from_arrays(geom, time, arrays))
    return out


# --------------------------------------------------------------------------- fake backend

_FAKES = weakref.WeakValueDictionary()
_tokens = itertools.count()


def _lookup_fake(token):
    try:
        return _FAKES[token]
    except KeyError:
        raise PhareseeIOError(f"fake backend {token} no longer exists") from None


class FakeBackend:
    """Backend protocol over in-memory data, for tests of lazy reading.

    data   : list of LevelData (one time) or ``{time: list of LevelData}``; the source values
    delay  : seconds to sleep inside every ``read_into`` (to make threads overlap)
    Counters: ``calls`` (number of read_into calls), ``requests`` (list of
    ``(time, level, file_ids)``), ``reads[(time, level)]`` (int64 count per file-order patch),
    ``index_calls``. ``fail_next(n)`` makes the next ``n`` reads write half of the patches then
    raise ``PhareseeIOError``. ``spec()`` only works inside this process (a registry token).
    """

    def __init__(self, data, quantity="fake", path="/nonexistent/fake.vtkhdf", delay=0.0):
        if not isinstance(data, dict):
            data = {data[0].time: data}
        self._data = {float(t): {ld.geom.level: ld for ld in lds} for t, lds in data.items()}
        first = next(iter(next(iter(self._data.values())).values()))
        g0 = first.geom
        r = 2**g0.level
        self._meta = FileMeta(
            path=path, format="fake", quantity=quantity, dim=g0.dim,
            dl0=np.asarray(g0.dl) * r, domain_cells0=np.asarray(g0.domain_cells) // r,
            periodic=np.asarray(g0.periodic), interp_order=1, model="hybrid",
            names=first.names, dtype=first.dtype, colocated=True, writer_signature="fake",
        )
        self.delay = float(delay)
        self.calls = 0
        self.index_calls = 0
        self.requests = []
        self.names_read = []  # component names of each read_into call (meta order)
        self.reads = {}
        self._fail = 0
        self._lock = threading.Lock()
        self.token = next(_tokens)
        _FAKES[self.token] = self

    def fail_next(self, n=1):
        self._fail = int(n)

    def meta(self):
        return self._meta

    def times(self):
        return np.array(sorted(self._data), dtype=np.float64)

    def _level(self, time, level):
        try:
            return self._data[float(time)][level]
        except KeyError:
            raise TimeNotFoundError(f"fake backend: no level {level} at time {time}") from None

    def index(self, time):
        self.index_calls += 1
        if float(time) not in self._data:
            raise TimeNotFoundError(f"fake backend: no time {time}")
        out = []
        for level in sorted(self._data[float(time)]):
            lo, hi = self._data[float(time)][level].geom.file_order_boxes()
            out.append(LevelIndex(level, lo, hi))
        return tuple(out)

    def nbytes(self, time, level, file_ids, names=None):
        ld = self._level(time, level)
        slots = ld.geom.store_index[np.asarray(file_ids, np.int64)]
        ncomp = len(ld.names) if names is None else len(tuple(names))
        return int(ld.geom.size[slots].sum()) * ld.dtype.itemsize * ncomp

    def read_into(self, time, level, file_ids, out, dst_start):
        src = self._level(time, level)
        file_ids = np.asarray(file_ids, np.int64)
        dst_start = np.asarray(dst_start, np.int64)
        if len(file_ids) != len(dst_start):
            raise ValueError("file_ids and dst_start differ in length")
        if len(file_ids) > 1 and (np.diff(file_ids) <= 0).any():
            raise ValueError("file_ids must be sorted and unique")
        if not out or not set(out) <= set(src.names):
            raise ValueError(f"out must hold components of {src.names}, got {tuple(out)}")
        names = tuple(n for n in src.names if n in out)
        with self._lock:
            self.calls += 1
            self.requests.append((float(time), level, file_ids.copy()))
            self.names_read.append(names)
            counts = self.reads.setdefault((float(time), level),
                                           np.zeros(src.geom.npatch, np.int64))
            counts[file_ids] += 1
            fail = self._fail > 0
            self._fail -= int(fail)
        if self.delay:
            _time.sleep(self.delay)
        geom = src.geom
        stop = len(file_ids) // 2 if fail else len(file_ids)
        for fid, dst in zip(file_ids[:stop], dst_start[:stop]):  # I/O scatter: per patch
            slot = geom.store_index[fid]
            a, b = geom.offset[slot], geom.offset[slot + 1]
            for name in names:
                out[name][dst : dst + (b - a)] = src._buf[name][a:b]
        if fail:
            raise PhareseeIOError(
                f"{self._meta.path}: simulated read failure at t={time}, level {level}, "
                f"{len(file_ids)} patches"
            )

    def spec(self):
        return (_lookup_fake, self.token)
