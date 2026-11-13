"""Reader of PHARE ``.vtkhdf`` files (VTK "OverlappingAMR", written by pharevtkhdf).

File structure (one file per quantity)::

    /                      attrs dimension, domain_box (upper cell, level 0),
                           periodicities (or boundary_conditions on older files),
                           interpOrder (hybrid only), pop_mass/pop_charge (populations)
    /VTKHDF                attrs Type="OverlappingAMR", Version, Origin
    /VTKHDF/Level<l>       attr Spacing (float32, dim)
        AMRBox             (nboxes, 6) int: lo0, up0, lo1, up1, lo2, up2 (inclusive cells)
        PointData/data     (nrows, ncomp): one row per node, components in columns
    /VTKHDF/Steps          attr NSteps; Values (time of each step)
        Level<l>/NumberOfAMRBox, AMRBoxOffset, PointDataOffset/data   (one entry per step)

The nodes of a patch (cell box [lower, upper], nodes lower..upper+1, x fastest)
are written ``copies`` times in a row (4 in 1D, 2 in 2D, 1 in 3D, so that VTK
sees 3D cells); the copies are identical and only the first one is read.

Reading strategies (DESIGN section 5.1), chosen per ``read_into`` call:

S1 "level span"  for large requests on files with small chunks (production
                 files: (1024, 1) chunks of 8 KiB): walk the HDF5 chunk B-tree to
                 the level above the leaves, read the byte span of the dump in
                 windows <= 256 MiB with parallel ``preadv``, locate the chunks by
                 parsing the leaves from memory, copy them out. Refused (falls
                 back to S3) on any file feature the walker does not know, and on
                 any inconsistency found while reading.
S2 "per patch"   small requests: ``dset[r0:r1, :]`` per patch, file order.
S3 "slab"        large requests otherwise: ``dset[a:b, :]`` in row windows <= 256 MiB.

Environment: ``PHARESEE_RAW_IO=0`` disables S1; ``PHARESEE_IO_THREADS`` sets the
number of reader threads (default 8).
"""

from __future__ import annotations

import bisect
import os
import threading
import warnings
from collections import OrderedDict
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from itertools import pairwise

import numpy as np

from ..errors import (
    FormatError,
    PhareseeFileChanged,
    PhareseeIOError,
    TimeNotFoundError,
)
from . import btree, rawio
from .h5open import open_h5
from .meta import FileMeta, LevelIndex

#: S1 and S3 read windows (bytes); S1 keeps two windows in memory (double buffering)
WINDOW_BYTES = 256 << 20
#: S1 only on files whose chunks are at most this many bytes
S1_MAX_CHUNK_BYTES = 64 << 10
#: a request is "large" (S1/S3) above this fraction of the level's bytes at that step
LARGE_FRACTION = 0.06
#: S3: two selected patches further apart than this are read by separate windows
S3_GAP_BYTES = 16 << 20
#: number of per-step index tables kept in memory
INDEX_CACHE_STEPS = 32

COPIES = {1: 4, 2: 2, 3: 1}
COMPONENT_NAMES = {
    1: ("value",),
    3: ("x", "y", "z"),
    6: ("xx", "xy", "xz", "yy", "yz", "zz"),
    9: ("xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy", "zz"),
}
_INT32_MAX = 2**31 - 1
#: S1 decode: chunk rows further apart than this are staged separately
_STAGE_GAP_CHUNKS = 64
_STRATEGIES = ("auto", "s1", "s2", "s3")


# ---------------------------------------------------------------- helpers


def round_shortest(values, dtype):
    """Round each value to the shortest decimal that round-trips `dtype`.

    float32 0.4 is stored as 0.4000000059604645; this returns 0.4 (float64).
    Returns a float64 array of the same shape.
    """
    dtype = np.dtype(dtype)
    arr = np.asarray(values, dtype=dtype)
    flat = [float(str(v)) for v in arr.reshape(-1)]
    return np.array(flat, dtype=np.float64).reshape(arr.shape)


def time_tolerance(t, dtype, given=None):
    """Matching tolerance at time `t`: ``rawio.time_tolerance`` (max(2 ulp of `dtype`
    at t, 5e-11), widened to 2 ulp of `given`'s dtype for a float32 request)."""
    return rawio.time_tolerance(t, dtype, given)


def _attr_str(value):
    if isinstance(value, bytes):
        return value.decode()
    if isinstance(value, np.ndarray) and value.shape == ():
        return _attr_str(value[()])
    return str(value)


def _exclusive_cumsum(a):
    out = np.zeros(len(a), dtype=np.int64)
    if len(a) > 1:
        np.cumsum(a[:-1], out=out[1:])
    return out


@dataclass(frozen=True)
class _CachedStep:
    """In-memory index of one step and what it was validated against.

    levels  tuple of LevelIndex (what ``index`` returns)
    sig     (st_size, st_mtime_ns) of the file when the entry was built or last re-validated
    prints  bytes: first and last data row of the step at every level (content fingerprint)
    """

    levels: tuple
    sig: tuple
    prints: bytes


@dataclass(frozen=True)
class _Payload:
    """Backend-private part of a LevelIndex.

    step       index of the step in the Steps tables
    row_start  (npatch,) int64 first row of each patch (first copy) in PointData/data
    npts       (npatch,) int64 number of nodes of each patch
    """

    step: int
    row_start: np.ndarray
    npts: np.ndarray


class _Torn(Exception):
    """The file changed while S1 was reading it."""


class _Scatter:
    """Copy contiguous file rows of one component into the requested patches.

    starts, npts : first-copy row range of each requested patch (sorted, disjoint)
    dst          : node offset of each patch in the output buffers
    outs         : one flat output buffer per component (file column order); None for a
                   column not asked for (``put`` is never called for it)
    One slice copy per (patch, put) intersection; ``count`` = elements written.
    """

    def __init__(self, starts, npts, dst, outs):
        self.starts = np.asarray(starts, dtype=np.int64)
        self.ends = self.starts + np.asarray(npts, dtype=np.int64)
        self._starts = self.starts.tolist()
        self._ends = self.ends.tolist()
        self._dst = np.asarray(dst, dtype=np.int64).tolist()
        self.outs = outs
        self.count = 0

    def put(self, comp, row_a, vals):
        """`vals` holds rows [row_a, row_a + len(vals)) of column `comp` (any stride)."""
        row_b = row_a + len(vals)
        starts, ends, dst = self._starts, self._ends, self._dst
        # bisect on lists: ~10x cheaper than np.searchsorted for one scalar
        p0 = bisect.bisect_right(ends, row_a)
        p1 = bisect.bisect_left(starts, row_b)
        out = self.outs[comp]
        for p in range(p0, p1):
            a = max(starts[p], row_a)
            b = min(ends[p], row_b)
            if b > a:
                d = dst[p] + a - starts[p]
                out[d : d + b - a] = vals[a - row_a : b - row_a]
                self.count += b - a


def _needed_chunk_rows(starts, npts, chunk_rows, k0, k1):
    """bool (k1-k0,): chunk rows that hold at least one requested row."""
    ka = starts // chunk_rows - k0
    kb = (starts + npts - 1) // chunk_rows - k0 + 1
    marks = np.zeros(k1 - k0 + 1, dtype=np.int64)
    np.add.at(marks, ka, 1)
    np.add.at(marks, kb, -1)
    return np.cumsum(marks[:-1]) > 0


def _windows(cuts, span_lo, span_hi, width):
    """Cut [span_lo, span_hi) into windows <= width, preferably at `cuts` (sorted)."""
    bounds = [span_lo]
    start = span_lo
    while start < span_hi:
        if span_hi - start <= width:
            stop = span_hi
        else:
            i = int(np.searchsorted(cuts, start + width, side="right")) - 1
            stop = int(cuts[i]) if i >= 0 and cuts[i] > start else start + width
        bounds.append(stop)
        start = stop
    return list(pairwise(bounds))


def _true_runs(mask):
    """(starts, stops) of the runs of True in a 1-D bool array."""
    edges = np.diff(np.concatenate([[0], mask.view(np.int8), [0]]))
    return np.nonzero(edges == 1)[0], np.nonzero(edges == -1)[0]


# ---------------------------------------------------------------- backend


class VtkhdfBackend:
    """Backend (see ``amr.io.Backend``) for one PHARE ``.vtkhdf`` quantity file.

    path          the .vtkhdf file
    strategy      "auto" (default) or "s1" / "s2" / "s3" to force a read strategy
                  ("s1" still falls back to S3 when the file does not allow it)
    window_bytes  S1/S3 window size in bytes (default WINDOW_BYTES = 256 MiB)
    identity      optional FileIdentity tuple recorded earlier (from ``spec()``):
                  raises PhareseeFileChanged if the file was replaced since

    Opening reads the attributes and the small Steps tables (a few ms), checks
    the format (FormatError naming the failed field) and the consistency of
    the Steps tables (DESIGN 5.1: an incomplete last dump is dropped with one
    RuntimeWarning). No file handle is kept: each call opens the file again.
    ``last_read`` describes the last read_into call (strategy, fallback reason).

    File changes: a file replaced (new inode) or shrunk raises
    PhareseeFileChanged. A file modified in place (same inode, other size or
    mtime) is accepted only if everything indexed so far is unchanged (a
    running simulation appending dumps): before a cached step index is used
    again, its Steps entries, its AMRBox rows and a fingerprint of its data
    (first and last row of the step at every level) are read again and
    compared. Any difference means the file was rewritten: PhareseeFileChanged,
    now and for every later call (a new backend reads the new content).
    """

    def __init__(self, path, *, strategy="auto", window_bytes=None, identity=None):
        if strategy not in _STRATEGIES:
            raise ValueError(f"strategy must be one of {_STRATEGIES}, not {strategy!r}")
        self._given_path = str(path)
        if identity is not None:
            identity = rawio.FileIdentity(*identity)
            rawio.check_identity(identity, path)
        self._identity = rawio.file_identity(path)
        self._opened = self._identity  # the Steps tables in memory describe this state
        self._path = self._identity.realpath
        self._strategy = strategy
        self._window = int(window_bytes) if window_bytes else WINDOW_BYTES
        self._lock = threading.Lock()
        self._index_cache = OrderedDict()
        self._rewritten = None  # PhareseeFileChanged once a rewrite was detected
        self.last_read = {}
        try:
            with self._open() as f:
                self._read_meta(f)
            self._s1_reasons = self._probe_s1()
        except (FormatError, PhareseeIOError):
            raise
        except (OSError, RuntimeError, KeyError) as e:
            # damaged metadata (e.g. a file truncated inside its headers) surfaces
            # as raw h5py errors while reading the attributes or tables
            raise PhareseeIOError(
                f"{self._path}: cannot read the file structure ({e}); the file is damaged "
                "or still being written"
            ) from e
        self._meta = self._build_meta()

    # ------------------------------------------------------------ opening

    def _open(self):
        return open_h5(self._path)

    def _fail(self, what):
        raise FormatError(f"{self._path}: {what}")

    def _read_meta(self, f):
        attrs = f.attrs
        if "VTKHDF" not in f:
            self._fail("no /VTKHDF group: not a vtkhdf file")
        root = f["VTKHDF"]
        kind = _attr_str(root.attrs.get("Type", ""))
        if kind != "OverlappingAMR":
            self._fail(f"/VTKHDF attribute Type is {kind!r}, expected 'OverlappingAMR'")
        if "dimension" not in attrs:
            self._fail("root attribute 'dimension' is missing")
        dim = int(np.asarray(attrs["dimension"]).reshape(-1)[0])
        if dim not in (1, 2, 3):
            self._fail(f"root attribute dimension = {dim}, expected 1, 2 or 3")
        if "domain_box" not in attrs:
            self._fail("root attribute 'domain_box' is missing")
        domain_box = np.asarray(attrs["domain_box"], dtype=np.int64).reshape(-1)
        if len(domain_box) != dim:
            self._fail(f"domain_box {domain_box.tolist()} does not have dimension={dim} values")
        self._dim = dim
        self._domain_cells = domain_box + 1
        self._periodic = self._read_periodicity(attrs, dim)
        self._interp = int(np.asarray(attrs["interpOrder"]).reshape(-1)[0]) if (
            "interpOrder" in attrs
        ) else None
        skip = {"dimension", "domain_box", "periodicities", "boundary_conditions", "interpOrder"}
        self._attrs = {k: self._plain(v) for k, v in attrs.items() if k not in skip}

        self._read_steps(f, root)

    @staticmethod
    def _plain(value):
        if isinstance(value, np.ndarray):
            return value.item() if value.size == 1 else value.tolist()
        if isinstance(value, np.generic):
            return value.item()
        if isinstance(value, bytes):
            return value.decode()
        return value

    def _read_periodicity(self, attrs, dim):
        if "periodicities" in attrs:
            per = np.asarray(attrs["periodicities"]).reshape(-1)
            if len(per) != dim:
                self._fail(f"periodicities {per.tolist()} does not have {dim} values")
            return per.astype(bool)
        if "boundary_conditions" in attrs:
            # DESIGN-NOTE: files written before the MHD boundary work (de97d2443) have
            # boundary_conditions = ["periodic", ...] instead of periodicities.
            bcs = [_attr_str(b) for b in np.asarray(attrs["boundary_conditions"]).reshape(-1)]
            if len(bcs) != dim:
                self._fail(f"boundary_conditions {bcs} does not have {dim} values")
            return np.array([b.strip().lower() == "periodic" for b in bcs], dtype=bool)
        self._fail("neither 'periodicities' nor 'boundary_conditions' root attribute found")

    def _read_steps(self, f, root):
        if "Steps" not in root or "Values" not in root["Steps"]:
            self._fail("no /VTKHDF/Steps/Values table")
        steps = root["Steps"]
        values_ds = steps["Values"]
        self._time_dtype = values_ds.dtype
        values = values_ds[()].reshape(-1)
        nsteps_attr = int(np.asarray(steps.attrs.get("NSteps", len(values))).reshape(-1)[0])

        levels = []
        while f"Level{len(levels)}" in root:
            levels.append(len(levels))
        if not levels:
            # DESIGN-NOTE: a file without data (momentum_tensor) fails at open, with a
            # FormatError naming what is missing, rather than at the first index().
            self._fail(
                "the file has no data (no /VTKHDF/Level0 group, only the Steps table): the "
                "writer created it but wrote no field (the vtkhdf writer does not write "
                "momentum_tensor); read this quantity from the phareh5 output instead"
            )
        self._levels = []
        tables = []
        dtype = ncomp = None
        for lvl in levels:
            info = self._read_level(root, steps, lvl)
            if dtype is None:
                dtype, ncomp = info["dtype"], info["ncomp"]
            elif info["dtype"] != dtype or info["ncomp"] != ncomp:
                self._fail(f"Level{lvl}/PointData/data differs from Level0 in dtype or columns")
            self._levels.append(info)
            tables += [info["nbox"], info["boxoff"], info["dataoff"]]
        self._dtype = np.dtype(dtype)
        self._names = COMPONENT_NAMES[ncomp]
        self._dl0 = round_shortest(self._levels[0]["spacing"], self._levels[0]["spacing_dtype"])

        nvalues = len(values)
        lengths = [len(t) for t in tables]
        m, big_m = min(lengths), max(lengths)
        if not (nvalues == nsteps_attr == m == big_m):
            # DESIGN-NOTE: NSteps only matters for the "normal" case; a level table
            # longer than Values (big_m > nvalues) is refused; the same warning is
            # given when only extra Values are ignored (emergency dump).
            tail_equal = nvalues > m and bool(np.all(values[m:] == values[m]))
            ok = big_m - m <= 1 and big_m <= nvalues and (nvalues - m <= 1 or tail_equal)
            if not ok:
                self._fail(
                    f"Steps tables disagree: Values has {nvalues} entries, NSteps={nsteps_attr}, "
                    f"level tables have {m} to {big_m} entries. Files written with per-level "
                    "dumps (fine_dump_lvl_max) are not supported"
                )
            if m < nvalues or big_m > m:
                warnings.warn(
                    f"{self._path}: last dump incomplete, ignored (Values has {nvalues} "
                    f"entries, level tables {m} to {big_m}); using the first {m} steps",
                    RuntimeWarning,
                    stacklevel=4,
                )
        self._nsteps = m
        self._raw_values = values[:m].copy()
        self._build_times(values[:m])

    def _read_level(self, root, steps, lvl):
        group = root[f"Level{lvl}"]
        where = f"/VTKHDF/Level{lvl}"
        if "PointData" not in group or "data" not in group["PointData"]:
            self._fail(f"{where}/PointData/data is missing")
        data = group["PointData/data"]
        if data.ndim != 2 or data.shape[1] not in COMPONENT_NAMES:
            self._fail(f"{where}/PointData/data has shape {data.shape}: expected (n, 1|3|6|9)")
        if data.dtype.kind != "f":
            self._fail(f"{where}/PointData/data has dtype {data.dtype}, expected floats")
        if "AMRBox" not in group:
            self._fail(f"{where}/AMRBox is missing")
        boxes = group["AMRBox"]
        # DESIGN-NOTE: 6 columns per the design; 2*dim columns is also accepted
        # (cheap, and the writer history mentions 2*dim-wide boxes).
        if boxes.ndim != 2 or boxes.shape[1] not in (6, 2 * self._dim):
            self._fail(f"{where}/AMRBox has shape {boxes.shape}: expected (n, 6)")
        if "Spacing" not in group.attrs:
            self._fail(f"{where} attribute Spacing is missing")
        spacing = np.asarray(group.attrs["Spacing"]).reshape(-1)
        if len(spacing) != self._dim:
            self._fail(f"{where} Spacing {spacing.tolist()} does not have {self._dim} values")
        if f"Level{lvl}" not in steps:
            self._fail(f"/VTKHDF/Steps/Level{lvl} is missing")
        st = steps[f"Level{lvl}"]
        tabs = {}
        for key in ("NumberOfAMRBox", "AMRBoxOffset", "PointDataOffset/data"):
            if key not in st:
                self._fail(f"/VTKHDF/Steps/Level{lvl}/{key} is missing")
            ds = st[key]
            if ds.dtype.kind != "i" or ds.dtype.itemsize not in (4, 8):
                self._fail(f"/VTKHDF/Steps/Level{lvl}/{key} has dtype {ds.dtype}: int32/int64")
            tabs[key] = (ds[()].reshape(-1).astype(np.int64), ds.dtype.itemsize == 4)
        return dict(
            level=lvl,
            dtype=data.dtype,
            ncomp=int(data.shape[1]),
            chunks=data.chunks,
            compression=data.compression,
            box_dtype=boxes.dtype,
            box_chunks=boxes.chunks,
            box_cols=int(boxes.shape[1]),
            spacing=spacing,
            spacing_dtype=spacing.dtype,
            nbox=tabs["NumberOfAMRBox"][0],
            boxoff=tabs["AMRBoxOffset"][0],
            dataoff=tabs["PointDataOffset/data"][0],
            int32=any(t[1] for t in tabs.values()),
            offset_dtype=st["AMRBoxOffset"].dtype,
        )

    def _build_times(self, raw):
        """Round, deduplicate (last copy wins), sort."""
        rounded = round_shortest(raw, self._time_dtype)
        # DESIGN-NOTE: the design dedupes consecutive equal times (restart); any
        # repeated time is handled the same way: the last copy wins (dict order).
        step_of = {t: s for s, t in enumerate(rounded.tolist())}
        times = np.array(sorted(step_of), dtype=np.float64)
        self._times = times
        self._steps = np.array([step_of[t] for t in times.tolist()], dtype=np.int64)

    def _probe_s1(self):
        """Per level: None if S1 can read it, else the reason it cannot."""
        reasons = []
        fd = rawio.open_readonly(self._path)
        try:
            with self._open() as f:
                for info in self._levels:
                    reasons.append(self._s1_refusal(f, fd, info))
        finally:
            os.close(fd)
        return reasons

    def _s1_refusal(self, f, fd, info):
        chunks = info["chunks"]
        if chunks is None:
            return "data not chunked"
        chunk_bytes = int(np.prod(chunks)) * info["dtype"].itemsize
        if chunk_bytes > S1_MAX_CHUNK_BYTES:
            return f"chunks of {chunk_bytes} bytes > {S1_MAX_CHUNK_BYTES} (h5py is fast enough)"
        try:
            btree.ChunkBTree(
                fd, self._header_addr(f, info["level"]), info["dtype"], chunks, info["ncomp"]
            )
        except btree.Unsupported as e:
            return str(e)
        return None

    @staticmethod
    def _header_addr(f, level):
        import h5py

        group = f[f"VTKHDF/Level{level}/PointData"]
        link = group.id.links.get_info(b"data")
        if link.type != h5py.h5l.TYPE_HARD:
            raise btree.Unsupported("PointData/data is not a hard link")
        return int(link.u)

    def _build_meta(self):
        lv0 = self._levels[0]
        s1 = self._s1_reasons[0]
        signature = (
            f"data {lv0['dtype']} chunks {lv0['chunks']}"
            + (f" {lv0['compression']}" if lv0["compression"] else "")
            + f", AMRBox {lv0['box_dtype']} x{lv0['box_cols']} chunks {lv0['box_chunks']}"
            + f", times {self._time_dtype}, offsets {lv0['offset_dtype']}"
            + (", raw reads (S1) on" if s1 is None else f", S1 off: {s1}")
        )
        name = os.path.basename(self._path)
        quantity = name[: -len(".vtkhdf")] if name.endswith(".vtkhdf") else name
        return FileMeta(
            path=self._path,
            format="vtkhdf",
            quantity=quantity,
            dim=self._dim,
            dl0=self._dl0,
            domain_cells0=self._domain_cells.copy(),
            periodic=self._periodic.copy(),
            interp_order=self._interp,
            model="hybrid" if self._interp is not None else "mhd",
            names=self._names,
            dtype=self._dtype,
            colocated=True,
            writer_signature=signature,
            attrs=dict(self._attrs),
        )

    # ------------------------------------------------------------ protocol

    def meta(self):
        """FileMeta of the file (read at open; no I/O)."""
        return self._meta

    def times(self):
        """Sorted unique dump times (float64, rounded to the file's time dtype). No I/O."""
        return self._times.copy()

    def spec(self):
        """("vtkhdf", path, kwargs): ``open_backend(path, **kwargs)`` reopens the file
        and raises PhareseeFileChanged if it was replaced in between."""
        return ("vtkhdf", self._path, {"identity": tuple(self._identity)})

    def _step_of(self, time):
        """Index into times() of `time` (tolerance: ``rawio.time_tolerance``, which
        also covers a float32 value passed by the caller)."""
        i = rawio.nearest_time(self._times, time, self._time_dtype)
        if i is not None:
            return i
        t = float(time)
        times = self._times
        order = np.argsort(np.abs(times - t), kind="stable")[:3]
        near = ", ".join(repr(float(times[j])) for j in sorted(order.tolist()))
        raise TimeNotFoundError(
            f"{self._path}: no dump at time {t!r}; nearest times: {near or 'none'}. "
            "List all times with times()."
        )

    def index(self, time):
        """Tuple of LevelIndex (FILE order) for the levels that have patches at `time`.

        Reads the AMRBox rows of that step only (one small read per level) and the
        first and last data row of the step per level (a content fingerprint); the
        result is cached. When the file was modified since the entry was made, the
        entry is checked against the file first (see the class docstring).
        Raises TimeNotFoundError, FormatError (inconsistent Steps offsets for that
        step), PhareseeFileChanged (file replaced, shrunk or rewritten in place).
        """
        return self._step_index(self._step_of(time))

    def _current_identity(self):
        """Stat the file: raises PhareseeFileChanged if it was replaced, shrunk or found
        rewritten earlier; returns the current FileIdentity."""
        if self._rewritten is not None:
            raise self._rewritten
        try:
            return rawio.check_identity(self._identity)
        except PhareseeFileChanged as e:  # replaced or shrunk: final as well
            with self._lock:
                self._rewritten = e
                self._index_cache.clear()
            raise

    def _mark_rewritten(self, what):
        """Forget every cached index and fail this and every later call."""
        err = rawio.rewritten(self._path, what)
        with self._lock:
            self._rewritten = err
            self._index_cache.clear()
        raise err

    def _step_index(self, i):
        now = self._current_identity()
        sig = (now.st_size, now.st_mtime_ns)
        with self._lock:
            hit = self._index_cache.get(i)
            if hit is not None:
                self._index_cache.move_to_end(i)
        if hit is not None and hit.sig == sig:
            return hit.levels
        step = int(self._steps[i])
        modified = rawio.modified_since(self._opened, now)
        with self._open() as f:
            if modified:
                self._check_steps(f, step, i)
            try:
                levels, prints = self._build_step(f, step)
            except FormatError:
                if modified:  # boxes no longer consistent with the Steps read at open
                    self._mark_rewritten(
                        f"the boxes of time {float(self._times[i])!r} no longer match"
                    )
                raise
        if hit is not None:
            if not self._same_step(hit, levels, prints):
                self._mark_rewritten(
                    f"the patches or the data of time {float(self._times[i])!r} differ"
                )
            levels = hit.levels  # unchanged: keep the objects already handed out
        with self._lock:
            if self._rewritten is not None:
                raise self._rewritten
            self._index_cache[i] = _CachedStep(levels, sig, prints)
            self._index_cache.move_to_end(i)
            while len(self._index_cache) > INDEX_CACHE_STEPS:
                self._index_cache.popitem(last=False)
            # the grown file is accepted (check_identity refused a shrink)
            self._identity = now
        return levels

    def _check_steps(self, f, step, i):
        """The Steps entries of `step` in the file still equal those read at open."""
        root = f["VTKHDF"]
        steps = root["Steps"]
        try:
            value = steps["Values"][step]
            same = value == self._raw_values[step]
            for info in self._levels:
                st = steps[f"Level{info['level']}"]
                for key, table in (("NumberOfAMRBox", "nbox"), ("AMRBoxOffset", "boxoff"),
                                   ("PointDataOffset/data", "dataoff")):
                    same = same and int(st[key][step]) == int(info[table][step])
        except (KeyError, ValueError, IndexError, OSError):
            same = False
        if not same:
            self._mark_rewritten(f"the Steps tables of time {float(self._times[i])!r} differ")

    def _build_step(self, f, step):
        """(tuple of LevelIndex, fingerprint bytes) of `step`, read from the open file."""
        out, prints = [], []
        for info in self._levels:
            lix = self._level_index(f, info, step)
            if lix is not None:
                out.append(lix)
                prints.append(self._fingerprint(f, lix))
        return tuple(out), b"".join(prints)

    @staticmethod
    def _fingerprint(f, lix):
        """Bytes of the first data row of the first patch and the last row of the last
        patch of a level at one step (two small reads). A pure append never changes
        them; a rewrite with other data almost surely does."""
        # DESIGN-NOTE: DESIGN 5 re-validates indexes on the Steps tables and boxes only;
        # a file rewritten with the same layout (same run setup, another quantity or
        # another seed) passes those, so two data rows of the step are compared too.
        data = f[f"VTKHDF/Level{lix.level}/PointData/data"]
        first = int(lix.payload.row_start[0])
        last = int(lix.payload.row_start[-1] + lix.payload.npts[-1] - 1)
        return data[first : first + 1].tobytes() + data[last : last + 1].tobytes()

    @staticmethod
    def _same_step(hit, levels, prints):
        if hit.prints != prints or len(hit.levels) != len(levels):
            return False
        for a, b in zip(hit.levels, levels):
            if (a.level != b.level or not np.array_equal(a.lower, b.lower)
                    or not np.array_equal(a.upper, b.upper)
                    or not np.array_equal(a.payload.row_start, b.payload.row_start)):
                return False
        return True

    def _level_index(self, f, info, s):
        lvl = info["level"]
        where = f"level {lvl}, step {s}"
        nbox, boxoff, dataoff = info["nbox"], info["boxoff"], info["dataoff"]
        nb, bo, po = int(nbox[s]), int(boxoff[s]), int(dataoff[s])
        if info["int32"]:
            # DESIGN-NOTE: saturation is checked on int32 tables only (an int64 table
            # cannot saturate), on entries s and s + 1 (those used for this step).
            nxt = slice(s, s + 2)
            if any((t[nxt] == _INT32_MAX).any() for t in (nbox, boxoff, dataoff)):
                self._fail(
                    f"{where}: a Steps offset equals 2**31-1 (int32 saturation of a file "
                    "written before 64-bit offsets): the data of this step cannot be located"
                )
        if nb < 0 or bo < 0 or po < 0:
            self._fail(f"{where}: negative Steps entry ({nb}, {bo}, {po})")
        if nb == 0:
            return None
        root = f["VTKHDF"]
        boxes_ds = root[f"Level{lvl}/AMRBox"]
        if bo + nb > boxes_ds.shape[0]:
            self._fail(
                f"{where}: AMRBox rows {bo}..{bo + nb} beyond the table ({boxes_ds.shape[0]})"
            )
        boxes = boxes_ds[bo : bo + nb, : 2 * self._dim].astype(np.int64)
        lower = np.ascontiguousarray(boxes[:, 0::2])
        upper = np.ascontiguousarray(boxes[:, 1::2])
        if (upper < lower).any():
            self._fail(f"{where}: AMRBox with upper < lower")
        npts = np.prod(upper - lower + 2, axis=1).astype(np.int64)
        rows = npts * COPIES[self._dim]
        total = int(rows.sum())
        if s + 1 < len(boxoff) and int(boxoff[s + 1]) != bo + nb:
            self._fail(f"{where}: AMRBoxOffset[{s + 1}] != AMRBoxOffset[{s}] + NumberOfAMRBox")
        if s + 1 < len(dataoff) and int(dataoff[s + 1]) != po + total:
            self._fail(
                f"{where}: PointDataOffset[{s + 1}] = {int(dataoff[s + 1])} != "
                f"PointDataOffset[{s}] + rows of the boxes = {po + total}"
            )
        nrows = root[f"Level{lvl}/PointData/data"].shape[0]
        if po + total > nrows:
            self._fail(f"{where}: rows {po}..{po + total} beyond PointData/data ({nrows} rows)")
        payload = _Payload(step=s, row_start=po + _exclusive_cumsum(rows), npts=npts)
        return LevelIndex(level=lvl, lower=lower, upper=upper, payload=payload)

    def _find_level(self, time, level):
        for lix in self.index(time):
            if lix.level == level:
                return lix
        raise ValueError(f"{self._path}: level {level} has no patch at time {time}")

    def nbytes(self, time, level, file_ids, names=None):
        """Bytes of the arrays read_into fills for these patches and components
        `names` (None: all; meta dtype). Uses the step index (read once and
        cached); no field bytes are read. Transient read windows (<= 2 x
        window_bytes) are not included."""
        lix = self._find_level(time, level)
        ids = np.asarray(file_ids, dtype=np.int64)
        ncomp = len(self._columns(names))
        return int(lix.payload.npts[ids].sum()) * ncomp * self._dtype.itemsize

    def _columns(self, names):
        """Sorted file columns of components `names` (None: all). Raises ValueError for
        an empty selection or an unknown name."""
        if names is None:
            return list(range(len(self._names)))
        names = [names] if isinstance(names, str) else list(names)
        unknown = [n for n in names if n not in self._names]
        if unknown or not names:
            raise ValueError(
                f"components {unknown or names} not among the file's {self._names}"
            )
        return sorted({self._names.index(n) for n in names})

    def reads_whole_level(self, time, level, file_ids):
        """True when ``read_into(time, level, file_ids, ...)`` would use S1 (level span).

        S1 reads every byte between the first and the last requested patch of
        the level, so the caller may as well ask for all the patches of that
        range (``LevelData.ensure`` then loads the whole level when it fits the
        memory budget). Uses the cached step index; no data is read.
        """
        lix = self._find_level(time, level)
        ids = np.asarray(file_ids, dtype=np.int64).reshape(-1)
        if len(ids) == 0:
            return False
        return self.choose_strategy(level, lix, lix.payload.npts[ids]) == "S1"

    @property
    def window_bytes(self):
        """Size of one S1/S3 read window; S1 holds two at a time."""
        return self._window

    def read_into(self, time, level, file_ids, out, dst_start):
        """Fill `out` with the patches `file_ids` of `level` at `time`.

        file_ids   sorted int64 FILE-order patch indices
        out        {component name: 1-D writable buffer} for a non-empty subset of
                   meta().names: only those columns are read (S2/S3 select them in
                   the file; S1 skips the chunks of the other columns)
        dst_start  int64 (len(file_ids),): node offset of each patch in the buffers

        Each patch is written x fastest, nodes lower..upper+1, in the buffers'
        dtype; only the first copy of each patch is read. Raises
        TimeNotFoundError, ValueError (bad arguments), FormatError,
        PhareseeFileChanged (file replaced/truncated), PhareseeIOError (read
        failure, naming file, time, level and patch count).
        """
        lix = self._find_level(time, level)
        ids = np.asarray(file_ids, dtype=np.int64).reshape(-1)
        dst = np.asarray(dst_start, dtype=np.int64).reshape(-1)
        if ids.shape != dst.shape:
            raise ValueError("read_into: file_ids and dst_start differ in length")
        npatch = len(lix.lower)
        if len(ids) and (ids[0] < 0 or ids[-1] >= npatch or (np.diff(ids) <= 0).any()):
            raise ValueError(f"read_into: file_ids must be sorted, unique, in [0, {npatch})")
        if not out or not set(out) <= set(self._names):
            raise ValueError(
                f"read_into: out must hold some of the components {self._names}, "
                f"got {tuple(out)}"
            )
        if len(ids) == 0:
            return
        starts = lix.payload.row_start[ids]
        npts = lix.payload.npts[ids]
        cols = self._columns(list(out))
        # one entry per file column, None for the columns not asked for
        outs = [out.get(n) for n in self._names]
        need = int((dst + npts).max())
        for c in cols:
            buf = outs[c]
            if buf.ndim != 1 or len(buf) < need:
                raise ValueError(
                    f"read_into: out[{self._names[c]!r}] must be 1-D with >= {need} elements"
                )
        rawio.check_identity(self._identity)
        try:
            self._dispatch(level, lix, starts, npts, dst, outs, cols)
        except (PhareseeFileChanged, FormatError):
            raise
        except (OSError, RuntimeError) as e:
            raise PhareseeIOError(
                f"{self._path}: reading time {float(time)!r}, level {level}, {len(ids)} "
                f"patches failed ({e}). Check that the file is readable and not being "
                "rewritten; PHARESEE_RAW_IO=0 disables the raw reader"
            ) from e

    # ------------------------------------------------------------ strategies

    def choose_strategy(self, level, lix, npts):
        """'S1', 'S2' or 'S3' for a request of `npts` nodes per patch (DESIGN 5.1)."""
        forced = self._strategy
        large = int(npts.sum()) > LARGE_FRACTION * int(lix.payload.npts.sum())
        raw_ok = rawio.raw_io_enabled()
        s1_ok = raw_ok and self._s1_reasons[level] is None
        if forced == "s2":
            return "S2"
        if forced == "s3":
            return "S3"
        if forced == "s1":
            return "S1" if raw_ok else "S3"
        if not large:
            return "S2"
        return "S1" if s1_ok else "S3"

    def _dispatch(self, level, lix, starts, npts, dst, outs, cols):
        strategy = self.choose_strategy(level, lix, npts)
        self.last_read = dict(strategy=strategy, fallback=None, npatch=len(starts),
                              columns=tuple(cols))
        with self._open() as f:
            dset = f[f"VTKHDF/Level{level}/PointData/data"]
            if strategy == "S2":
                return self._read_s2(dset, starts, npts, dst, outs, cols)
            if strategy == "S1":
                reason = self._read_s1(f, dset, level, starts, npts, dst, outs, cols)
                if reason is None:
                    return None
                self.last_read.update(strategy="S3", fallback=reason)
            return self._read_s3(dset, starts, npts, dst, outs, cols)

    def _read_s2(self, dset, starts, npts, dst, outs, cols):
        """Per patch: dset[r0:r1, :] in file order (dset[r0:r1, c] per wanted column c
        when not all columns are wanted).

        Plain slicing, not ``read_direct``: h5py serves simple slices through its
        compiled fast reader, ~3x faster per patch than ``read_direct`` on
        run0055Lx (143 patches: 1.1 ms for one column, 3.2 ms for three, against
        8.1 / 10.4 ms with ``read_direct``)."""
        if len(cols) < dset.shape[1]:
            for s, n, d in zip(starts.tolist(), npts.tolist(), dst.tolist()):
                for c in cols:
                    outs[c][d : d + n] = dset[s : s + n, c]
            return
        for s, n, d in zip(starts.tolist(), npts.tolist(), dst.tolist()):
            buf = dset[s : s + n, :]
            for c, out in enumerate(outs):
                out[d : d + n] = buf[:, c]

    def _read_s3(self, dset, starts, npts, dst, outs, cols):
        """Row windows <= window_bytes covering the requested patches, then scatter
        (one column selection per wanted column when not all columns are wanted)."""
        ncol = dset.shape[1]
        row_bytes = ncol * dset.dtype.itemsize
        wrows = max(1, self._window // row_bytes)
        gap_rows = max(1, S3_GAP_BYTES // row_bytes)
        total_rows = int((starts[-1] + npts[-1]) - starts[0])
        subset = len(cols) < ncol
        # (rows, 1) buffer for one column: same shape as the file selection (fast path)
        scratch = np.empty((min(wrows, total_rows), 1 if subset else ncol), dtype=dset.dtype)
        scatter = _Scatter(starts, npts, dst, outs)

        def flush(a, b):
            for x in range(a, b, wrows):
                y = min(b, x + wrows)
                buf = scratch[: y - x]
                if subset:
                    for c in cols:
                        dset.read_direct(buf, np.s_[x:y, c : c + 1])
                        scatter.put(c, x, buf[:, 0])
                    continue
                dset.read_direct(buf, np.s_[x:y, :])
                for c in range(ncol):
                    scatter.put(c, x, buf[:, c])

        cur_a = cur_b = None
        for s, n in zip(starts.tolist(), npts.tolist()):
            if cur_a is None:
                cur_a, cur_b = s, s + n
            elif s - cur_b > gap_rows or s + n - cur_a > wrows:
                flush(cur_a, cur_b)
                cur_a, cur_b = s, s + n
            else:
                cur_b = s + n
        flush(cur_a, cur_b)
        if scatter.count != len(cols) * int(npts.sum()):
            raise RuntimeError("S3 wrote an unexpected number of values")

    # S1 ------------------------------------------------------------------

    def _read_s1(self, f, dset, level, starts, npts, dst, outs, cols):
        """Try S1 (twice if the file changes during the read). Returns None on
        success, else the reason to fall back to S3."""
        reason = None
        for _attempt in range(2):
            try:
                self._s1_once(f, dset, level, starts, npts, dst, outs, cols)
                return None
            except _Torn:
                reason = "file modified during the raw read"
            except btree.Unsupported as e:
                return f"S1 refused: {e}"
        return reason

    def _s1_once(self, f, dset, level, starts, npts, dst, outs, cols):
        ncol = int(dset.shape[1])
        dtype = dset.dtype.newbyteorder("<")
        header = self._header_addr(f, level)
        fd = rawio.open_readonly(self._path)
        try:
            before = rawio.stat_signature(fd)
            bt = btree.ChunkBTree(fd, header, dset.dtype, dset.chunks, ncol)
            reader = _S1Reader(self, fd, bt, dtype, starts, npts, dst, outs, cols)
            nthreads = rawio.io_threads()
            with ThreadPoolExecutor(max_workers=nthreads) as pool, ThreadPoolExecutor(
                max_workers=1
            ) as prefetch:
                reader.run(pool, prefetch, self._window)
            if rawio.stat_signature(fd) != before:
                raise _Torn()
        finally:
            os.close(fd)
        self.last_read.update(reader.stats)


class _S1Reader:
    """One S1 read of rows of one level dataset (see module docstring)."""

    def __init__(self, backend, fd, bt, dtype, starts, npts, dst, outs, cols):
        self.fd, self.bt, self.dtype = fd, bt, dtype
        self.ncol = bt.ncols
        # chunks of the columns not asked for are dropped from the chunk tables (never
        # copied; read from disk only when they lie inside a window read anyway)
        self.cols = np.asarray(cols, dtype=np.int64)
        self.colmask = np.zeros(self.ncol, dtype=bool)
        self.colmask[self.cols] = True
        self.subset = len(self.cols) < self.ncol
        self.ch = bt.chunk_rows
        self.cb = bt.chunk_bytes
        self.scatter = _Scatter(starts, npts, dst, outs)
        self.npts_total = int(npts.sum())
        row_lo = int(starts[0])
        row_hi = int(starts[-1] + npts[-1])
        self.row_lo, self.row_hi = row_lo, row_hi
        self.k0 = row_lo // self.ch
        self.k1 = (row_hi - 1) // self.ch + 1
        self.needed = _needed_chunk_rows(starts, npts, self.ch, self.k0, self.k1)
        self.decoded = []  # linear chunk ids (krel * ncol + col) decoded
        self._stage = np.empty(0, dtype=dtype)  # reused across windows (no new page faults)
        self.stats = dict(windows=0, strays=0, leaves=0, bytes_read=0, t_read=0.0, t_decode=0.0)

    # chunk tables are (k, c, addr) int64 triples, k = absolute chunk row
    def _keep(self, rows, cols, addrs):
        k = rows // self.ch
        inside = (k >= self.k0) & (k < self.k1)
        if self.subset:
            inside &= self.colmask[cols]
        k, cols, addrs = k[inside], cols[inside], addrs[inside]
        sel = self.needed[k - self.k0]
        return k[sel], cols[sel], addrs[sel]

    def run(self, pool, prefetch, window):
        import time

        bt, fd = self.bt, self.fd
        span = bt.level_span(fd, self.row_lo, self.row_hi, executor=pool)
        main = span["main"]
        self.stats["leaves"] = len(span["leaves"])
        pend = self._keep(*bt.read_leaves(fd, span["others"], executor=pool))
        windows = _windows(main, span["span_lo"], span["span_hi"], max(window, 4 * bt.nodesize))
        parsed = np.zeros(len(main), dtype=bool)
        strays = []
        width = max(b - a for a, b in windows)
        bufs = [np.empty(width, dtype=np.uint8), np.empty(width, dtype=np.uint8)]

        def read(i):
            a, b = windows[i]
            buf = bufs[i % 2][: b - a]
            rawio.read_ranges(fd, [a], [b - a], out=buf, executor=pool)
            return buf

        fut = prefetch.submit(read, 0)
        for i, (w0, w1) in enumerate(windows):
            t0 = time.perf_counter()
            buf = fut.result()
            self.stats["t_read"] += time.perf_counter() - t0
            self.stats["bytes_read"] += w1 - w0
            if i + 1 < len(windows):
                fut = prefetch.submit(read, i + 1)
            t0 = time.perf_counter()
            inwin = (main >= w0) & (main + bt.nodesize <= w1)
            parsed |= inwin
            k, c, a = self._keep(*bt.parse_leaves(buf, main[inwin] - w0))
            k = np.concatenate([pend[0], k])
            c = np.concatenate([pend[1], c])
            a = np.concatenate([pend[2], a])
            here = (a >= w0) & (a + self.cb <= w1)
            later = a >= w1
            gone = ~(here | later)
            strays.append((k[gone], c[gone], a[gone]))
            pend = (k[later], c[later], a[later])
            self._decode(buf, k[here], c[here], a[here] - w0)
            self.stats["windows"] += 1
            self.stats["t_decode"] += time.perf_counter() - t0
        # leaves straddling a window cut, chunks behind us or beyond the span
        late = self._keep(*bt.read_leaves(fd, main[~parsed], executor=pool))
        k = np.concatenate([pend[0], late[0]] + [s[0] for s in strays])
        c = np.concatenate([pend[1], late[1]] + [s[1] for s in strays])
        a = np.concatenate([pend[2], late[2]] + [s[2] for s in strays])
        self.stats["strays"] = len(a)
        if len(a):
            t0 = time.perf_counter()
            order = np.argsort(a, kind="stable")
            k, c, a = k[order], c[order], a[order]
            buf = rawio.read_ranges(fd, a, np.full(len(a), self.cb), executor=pool)
            self.stats["bytes_read"] += len(a) * self.cb
            self._decode(buf, k, c, np.arange(len(a), dtype=np.int64) * self.cb)
            self.stats["t_read"] += time.perf_counter() - t0
        self._verify()

    def _decode(self, buf, k, c, rel):
        """Copy chunks (chunk row k, column c, at byte rel of buf) into the outputs.

        Contiguous runs of chunks (consecutive in the buffer and in (row, col)
        order) are copied per component with one strided copy into a row
        staging array, then each contiguous staged row range goes to the
        patches it covers (one slice copy per patch).
        """
        if len(k) == 0:
            return
        # chunks far apart in rows (strays) are staged separately: the staging
        # array spans the rows of one cluster only
        by_row = np.argsort(k, kind="stable")
        cuts = np.nonzero(np.diff(k[by_row]) > _STAGE_GAP_CHUNKS)[0] + 1
        if len(cuts):
            for part in np.split(by_row, cuts):
                self._decode_cluster(buf, k[part], c[part], rel[part])
        else:
            self._decode_cluster(buf, k, c, rel)

    def _decode_cluster(self, buf, k, c, rel):
        if self.subset:
            return self._decode_cluster_cols(buf, k, c, rel)
        ncol, ch, cb = self.ncol, self.ch, self.cb
        order = np.argsort(rel, kind="stable")
        k, c, rel = k[order], c[order], rel[order]
        j = (k - self.k0) * ncol + c
        self.decoded.append(j)
        brk = np.nonzero((np.diff(rel) != cb) | (np.diff(j) != 1))[0] + 1
        run_a = np.concatenate([[0], brk]).tolist()
        run_b = np.concatenate([brk, [len(k)]]).tolist()
        kmin = int(k.min())
        nk = int(k.max()) - kmin + 1
        if self._stage.size < ncol * nk * ch:
            self._stage = np.empty(ncol * nk * ch, dtype=self.dtype)
        stage = self._stage[: ncol * nk * ch].reshape(ncol, nk, ch)
        filled = np.zeros((ncol, nk), dtype=bool)
        jl, rl = j.tolist(), rel.tolist()
        base = kmin - self.k0
        for a, b in zip(run_a, run_b):
            m, j0 = b - a, jl[a]
            block = buf[rl[a] : rl[a] + m * cb].view(self.dtype).reshape(m, ch)
            for comp in range(ncol):
                i0 = (comp - j0) % ncol
                if i0 >= m:
                    continue
                part = block[i0::ncol]
                kk = (j0 + i0) // ncol - base
                stage[comp, kk : kk + len(part)] = part
                filled[comp, kk : kk + len(part)] = True
        for comp in range(ncol):
            ra, rb = _true_runs(filled[comp])
            for x, y in zip(ra.tolist(), rb.tolist()):
                self.scatter.put(comp, (kmin + x) * ch, stage[comp, x:y].reshape(-1))

    def _decode_cluster_cols(self, buf, k, c, rel):
        """``_decode_cluster`` when only some columns are wanted: per column, runs of
        chunks of consecutive rows lying at a constant byte stride in the buffer (e.g.
        ncol chunks apart when the columns are interleaved) are copied with one strided
        view each."""
        ncol, ch, cb = self.ncol, self.ch, self.cb
        itemsize = self.dtype.itemsize
        self.decoded.append((k - self.k0) * ncol + c)
        for comp in np.unique(c).tolist():
            sel = c == comp
            kc, rc = k[sel], rel[sel]
            order = np.argsort(kc, kind="stable")
            kc, rc = kc[order], rc[order]
            dk, dr = np.diff(kc), np.diff(rc)
            step = dk == 1
            stride = cb
            if step.any():
                vals, counts = np.unique(dr[step], return_counts=True)
                stride = int(vals[np.argmax(counts)])
            if stride < cb:  # chunks overlapping in the buffer: one chunk per run
                brk = np.arange(1, len(kc))
                stride = cb
            else:
                brk = np.nonzero(~step | (dr != stride))[0] + 1
            run_a = np.concatenate([[0], brk]).tolist()
            run_b = np.concatenate([brk, [len(kc)]]).tolist()
            kl, rl = kc.tolist(), rc.tolist()
            for a, b in zip(run_a, run_b):
                block = np.ndarray((b - a, ch), dtype=self.dtype, buffer=buf, offset=rl[a],
                                   strides=(stride, itemsize))
                self.scatter.put(comp, kl[a] * ch, block.reshape(-1))

    def _verify(self):
        """Every needed chunk decoded exactly once, every value written once."""
        got = np.sort(np.concatenate(self.decoded)) if self.decoded else np.empty(0, np.int64)
        krel = np.nonzero(self.needed)[0]
        expect = (krel[:, None] * self.ncol + self.cols).reshape(-1)
        if len(got) != len(expect) or not np.array_equal(got, expect):
            raise btree.Unsupported(
                f"chunk index inconsistent with the request ({len(got)} chunks found, "
                f"{len(expect)} expected)"
            )
        if self.scatter.count != len(self.cols) * self.npts_total:
            raise btree.Unsupported("raw read wrote an unexpected number of values")
