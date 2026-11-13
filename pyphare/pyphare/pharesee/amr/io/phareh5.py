"""Reader of PHARE's ``phareh5`` field files (one quantity per file).

File layout (written by src/diagnostic/detail/h5writer.hpp)::

    /                      attrs cell_width, dimension, domain_box, interpOrder
                           (hybrid only), layoutType, origin, periodicities,
                           (pop files) pop_mass, pop_charge
    /t/<"%.10f" time>/pl<level>/p<rank>#<localid>/      attrs lower, upper, ...
                                                 <dataset>  attr ghosts

Each dataset is one component of one patch, contiguous, C order ``[x, y, z]``,
with ``ghosts`` ghost nodes on each side, Yee-centred for hybrid quantities
(cell-centred for MHD fluid quantities).

What this module gives (DESIGN section 5.2):
- ``H5Backend.index(time)``: per level, the patch boxes in FILE order (patch
  groups sorted numerically by ``(rank, localid)``) and per dataset its ghost
  width, shape, byte offset and dtype. Built by a walk with the low-level h5py
  API, cached in memory and on disk (``indexcache``).
- ``H5Backend.read_into(...)``: reads datasets with threaded ``os.preadv`` at
  their byte offsets (``read_direct`` when a dataset has no single offset),
  converts each patch to node-centred values without ghosts (the same 2-point
  averages as the C++ vtkhdf writer) and writes them x fastest into flat
  buffers.
- ``H5Backend.read_derived("J_yee" | "divB_yee", ...)``: master's Yee-stencil
  current and divergence, computed per patch from the raw Yee B, averaged to
  nodes.

h5py is imported inside functions only (it may pull a different libhdf5).
"""

from __future__ import annotations

import os
import threading
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass

import numpy as np

from ..errors import (
    FormatError,
    PhareseeFileChanged,
    PhareseeIOError,
    TimeNotFoundError,
)
from . import FileMeta, LevelIndex, indexcache, rawio
from .h5open import open_h5
from .rawio import check_identity, file_identity

# counters of this process; tests read STATS["walks"] (time-group walks done)
# to check that a cached index is used instead of a walk
STATS = {"walks": 0}

BATCH_BYTES = 64 << 20  # raw bytes read per batch (two batches in flight)
COALESCE_GAP = 64 << 10  # datasets closer than this are read with one pread
PREFETCH_WINDOW = 256 << 20  # scratch window of the page-cache warm-up
PREFETCH_FRACTION = 0.5  # read_into prefetches when it asks for >= this share of a time


# dataset name -> key of pyphare.core.gridlayout.yee_centering.
# Same mapping as pharesee/hierarchy/hierarchy_utils.py:17-52 (field_qties),
# copied because importing that module imports the whole legacy hierarchy.
_HYBRID_KEYS = {
    "EM_B_x": "Bx",
    "EM_B_y": "By",
    "EM_B_z": "Bz",
    "EM_E_x": "Ex",
    "EM_E_y": "Ey",
    "EM_E_z": "Ez",
    "flux_x": "Fx",
    "flux_y": "Fy",
    "flux_z": "Fz",
    "bulkVelocity_x": "Vx",
    "bulkVelocity_y": "Vy",
    "bulkVelocity_z": "Vz",
    "charge_flux_x": "CFx",
    "charge_flux_y": "CFy",
    "charge_flux_z": "CFz",
    "momentum_tensor_xx": "Mxx",
    "momentum_tensor_xy": "Mxy",
    "momentum_tensor_xz": "Mxz",
    "momentum_tensor_yy": "Myy",
    "momentum_tensor_yz": "Myz",
    "momentum_tensor_zz": "Mzz",
    "density": "rho",
    "mass_density": "rho",
    "charge_density": "rho",
}
# MHD fluid quantities are cell-centred (gridlayout_mhd_yee.hpp:60-76);
# B is face-centred as in hybrid
_MHD_KEYS = {
    "EM_B_x": "Bx",
    "EM_B_y": "By",
    "EM_B_z": "Bz",
    "EM_E_x": "Ex",
    "EM_E_y": "Ey",
    "EM_E_z": "Ez",
    "rho": "mhdRho",
    "V_x": "mhdVx",
    "V_y": "mhdVy",
    "V_z": "mhdVz",
    "P": "mhdP",
    "rhoV_x": "mhdRhoVx",
    "rhoV_y": "mhdRhoVy",
    "rhoV_z": "mhdRhoVz",
    "Etot": "mhdEtot",
}
_VECTOR = ("x", "y", "z")
_TENSOR = ("xx", "xy", "xz", "yy", "yz", "zz")
_B_DATASETS = ("EM_B_x", "EM_B_y", "EM_B_z")


def dual_axes(model: str, dataset_name: str, dim: int) -> tuple:
    """Tuple of `dim` bools: True where the dataset is dual (cell-centred) along x, y, z.

    `model` is "hybrid" or "mhd". Uses pyphare.core.gridlayout.yee_centering.
    Raises FormatError for a dataset name this reader does not know.
    """
    from pyphare.core.gridlayout import yee_centering

    keys = _MHD_KEYS if model == "mhd" else _HYBRID_KEYS
    if dataset_name not in keys:
        raise FormatError(
            f"dataset {dataset_name!r} ({model} file): unknown centering; "
            f"known names are {sorted(keys)}. Use Run(path, legacy=True) for this file."
        )
    key = keys[dataset_name]
    return tuple(yee_centering[d][key] == "dual" for d in "xyz"[:dim])


def component_names(dataset_names) -> tuple:
    """(component names, dataset names in the same order) for the datasets of one patch.

    Vector quantities (datasets ending in _x, _y, _z) give ("x", "y", "z"), the
    momentum tensor gives ("xx", "xy", "xz", "yy", "yz", "zz"), a single
    dataset gives ("value",). Raises FormatError otherwise.
    """
    dataset_names = list(dataset_names)
    if len(dataset_names) == 1:
        return ("value",), (dataset_names[0],)
    for names in (_VECTOR, _TENSOR):
        if len(dataset_names) != len(names):
            continue
        found = {}
        for d in dataset_names:
            stem, _, suffix = d.rpartition("_")
            if suffix in names and stem:
                found[suffix] = (stem, d)
        stems = {stem for stem, _ in found.values()}
        if len(found) == len(names) and len(stems) == 1:
            return names, tuple(found[n][1] for n in names)
    raise FormatError(
        f"datasets {dataset_names}: expected one scalar, a vector (_x, _y, _z) or the "
        "momentum tensor (_xx ... _zz)"
    )


def plan_coalesced(offsets, nbytes, gap=COALESCE_GAP):
    """Group byte ranges into contiguous file runs.

    offsets, nbytes : int64 arrays (any order), ranges to read
    Returns (starts, stops, buf_starts, pos, total): file runs [start, stop) in
    increasing offset order, where each run lands in the staging buffer
    (buf_starts), the staging position of each input range (pos, same order as
    the input) and the staging size. Ranges separated by at most `gap` bytes
    share one run (the gap bytes are read and ignored).
    """
    offsets = np.asarray(offsets, np.int64)
    nbytes = np.asarray(nbytes, np.int64)
    if offsets.size == 0:
        e = np.zeros(0, np.int64)
        return e, e, e, e, 0
    order = np.argsort(offsets, kind="stable")
    o = offsets[order]
    end = o + nbytes[order]
    run_end = np.maximum.accumulate(end)
    new_run = np.ones(len(o), bool)
    new_run[1:] = o[1:] > run_end[:-1] + gap
    run_id = np.cumsum(new_run) - 1
    first = np.flatnonzero(new_run)
    starts = o[first]
    last = np.append(first[1:], len(o)) - 1
    stops = run_end[last]
    sizes = stops - starts
    buf_starts = np.cumsum(sizes) - sizes
    pos = np.empty(len(o), np.int64)
    pos[order] = buf_starts[run_id] + (o - starts[run_id])
    return starts, stops, buf_starts, pos, int(sizes.sum())


def _crop(arr, g, nn, dual):
    """Sub-array of the raw patch array needed for nodes 0..nn-1: one ghost cell
    on each dual side (the layer next to the domain), no ghost on primal axes."""
    sl = tuple(
        slice(g - 1, g + n) if d else slice(g, g + n) for n, d in zip(nn, dual)
    )
    return arr[sl]


def yee_to_nodes(arr, g, nn, dual, dst):
    """Average one raw Yee patch array to its nodes and write it into `dst`.

    arr  : raw dataset, indexed [x, y, z] (C order as in the file), with `g`
           ghost nodes per side
    g    : ghost width; at least 2 where `dual` holds (asserted by the caller):
           only the ghost layer next to the domain is used, never the outer one
    nn   : node counts per axis (upper - lower + 2)
    dual : bools per axis
    dst  : array of shape `nn` (usually a Fortran-order view into a flat
           buffer, i.e. x fastest)
    Along each dual axis node i = 0.5 * (v[i-1] + v[i]) (cells around the node),
    nested over dual axes in x, y, z order. One dual axis: computed in the
    file's dtype (one rounding, the correctly rounded mean). Two or more dual
    axes: intermediate averages in float64, one rounding at the end.
    """
    # DESIGN-NOTE: DESIGN 5.2 says "arithmetic in the file's dtype". With 2+
    # dual axes (Bz in 2D, B in 3D, MHD cell quantities) nested float32
    # averages round twice and lose up to ~100 ulp on small values next to
    # large ones (cancellation; 4.4 % of Bz nodes of sample_data/2d_static
    # exceed 1 ulp vs master's float64). Intermediates in float64 give the
    # correctly rounded value: <= 1 ulp from master, as the C++ vtkhdf writer
    # (which averages in double).
    x = _crop(arr, g, nn, dual)
    axes = [a for a, d in enumerate(dual) if d]
    if not axes:
        dst[...] = x
        return
    if len(axes) > 1 and x.dtype.itemsize < 8:
        x = x.astype(np.float64)
    # sums first, one scaling by 0.5**len(axes) at the end: a power of two, so
    # the result is the same as halving after each sum, with fewer passes
    for i, a in enumerate(axes):
        lo = x[(slice(None),) * a + (slice(0, -1),)]
        hi = x[(slice(None),) * a + (slice(1, None),)]
        if i == len(axes) - 1:
            np.add(lo, hi, out=dst)
        else:
            x = lo + hi
    np.multiply(dst, 0.5 ** len(axes), out=dst)


def yee_current(B, g, ncells, dl):
    """Master's Yee current (pharesee/run/utils.py:6-115), any dimension, one patch.

    B      : (Bx, By, Bz) raw Yee arrays [x, y, z] with `g` ghosts
    ncells : cells per axis; dl : float64 mesh size per axis
    Returns (Jx, Jy, Jz) float64 at their Yee positions with ghosts (J_a dual
    along a, primal elsewhere). Differences are taken in the file's dtype then
    divided by dl in float64, as master does. The outermost layer of each
    array along a derivative axis is NaN (master copies its neighbour; it is
    never used for nodes).
    """
    dim = len(ncells)
    out = []
    for a in range(3):
        shape = tuple(int(n) + 2 * g + (0 if d == a else 1) for d, n in enumerate(ncells))
        b, c = (a + 1) % 3, (a + 2) % 3
        terms = [(1.0, b, c)] if b < dim else []
        if c < dim:
            terms.append((-1.0, c, b))
        if not terms:
            out.append(np.zeros(shape))
            continue
        J = np.full(shape, np.nan)
        deriv_axes = {t[1] for t in terms}
        region = tuple(slice(1, -1) if d in deriv_axes else slice(None) for d in range(dim))
        acc = None
        for sign, axis, comp in terms:
            src = B[comp]
            hi = src[(slice(None),) * axis + (slice(1, None),)]
            lo = src[(slice(None),) * axis + (slice(0, -1),)]
            term = (hi - lo) / dl[axis]
            idx = tuple(
                slice(None) if (d == axis or d not in deriv_axes) else slice(1, -1)
                for d in range(dim)
            )
            term = term[idx]
            if sign < 0:
                term = -term
            acc = term if acc is None else acc + term
        J[region] = acc
        out.append(J)
    return tuple(out)


def yee_divergence(B, dl):
    """Master's Yee divB (pharesee/run/utils.py:118-121), any dimension, one patch.

    Returns float64 cell values (dual along every axis), with ghosts:
    sum over axes a < dim of (B_a[i+1] - B_a[i]) / dl[a], differences in the
    file's dtype.
    """
    dim = len(dl)
    acc = None
    for a in range(dim):
        src = B[a]
        hi = src[(slice(None),) * a + (slice(1, None),)]
        lo = src[(slice(None),) * a + (slice(0, -1),)]
        term = (hi - lo) / dl[a]
        acc = term if acc is None else acc + term
    return acc


def _parse_patch_name(name: bytes):
    text = name.decode() if isinstance(name, bytes) else str(name)
    try:
        rank, _, local = text[1:].partition("#")
        if not text.startswith("p") or not local:
            raise ValueError
        return int(rank), int(local)
    except ValueError:
        raise FormatError(
            f"patch group {text!r}: expected a name p<rank>#<localid>"
        ) from None


def _parse_level_name(name: bytes) -> int:
    text = name.decode() if isinstance(name, bytes) else str(name)
    if not text.startswith("pl") or not text[2:].isdigit():
        raise FormatError(f"level group {text!r}: expected a name pl<level>")
    return int(text[2:])


@dataclass(frozen=True)
class H5LevelTable:
    """Index of one level at one time (LevelIndex.payload of H5Backend), FILE order.

    rank, localid : (npatch,) int64 parsed from the patch group name
    ghosts        : (ncomp, npatch) int64 ghost width of each dataset
    shape         : (ncomp, npatch, dim) int64 dataset shape [x, y, z]
    offset        : (ncomp, npatch) int64 byte offset in the file, -1 if undefined
    nbytes        : (ncomp, npatch) int64 dataset size in bytes
    datasets      : dataset name of each component
    dtypes        : numpy dtype of each component in the file
    """

    time_key: str
    level: int
    rank: np.ndarray
    localid: np.ndarray
    ghosts: np.ndarray
    shape: np.ndarray
    offset: np.ndarray
    nbytes: np.ndarray
    datasets: tuple
    dtypes: tuple

    def patch_name(self, i) -> str:
        """HDF5 group name of file-order patch i."""
        return f"p{int(self.rank[i])}#{int(self.localid[i])}"


@dataclass(frozen=True)
class _CachedTime:
    """In-memory index of one time and what it was validated against.

    levels  {level: LevelIndex}
    sig     (st_size, st_mtime_ns) of the file when built or last re-validated
    digest  sha256 of the bytes at the t/<time> object header (at most 4 KiB, never
            past the next time's header as seen at open)
    prints  bytes: first 64 bytes of the first and of the last dataset (by offset) of
            every level (content fingerprint; empty for datasets without offset)
    """

    levels: dict
    sig: tuple
    digest: str
    prints: bytes


_PRINT_BYTES = 64


def _open_h5(path):
    """Open `path` read-only with h5py (``h5open.open_h5``: no file locking, or h5py's
    default locking when another handle of this process holds the file that way);
    PhareseeIOError naming the file if that fails."""
    return open_h5(path)


class H5Backend:
    """Backend (DESIGN section 5) for one phareh5 file.

    path      : the .h5 file of one quantity
    nthreads  : threads for preadv (default env PHARESEE_IO_THREADS or 8)
    raw_io    : False reads every dataset with h5py read_direct (also env
                PHARESEE_RAW_IO=0); results are identical, only slower
    prefetch  : "auto" (read_into asking for >= 50 % of a time's bytes warms the
                page cache over that time's byte region first), "never", or
                "walk" (also before an index walk that the cache cannot avoid;
                good for cold whole-domain work on slow file systems)
    Opening reads the root attributes, the list of times and the dataset names
    of one patch (one open of the file, closed before returning). Raises
    FormatError when the file is not a phareh5 field file, PhareseeIOError
    (naming the file) when it cannot be read as HDF5 (empty, truncated),
    ValueError for a bad ``PHARESEE_IO_THREADS`` / ``PHARESEE_RAW_IO`` value.

    File changes: a file replaced (new inode) or shrunk raises
    PhareseeFileChanged. A file modified in place (same inode, other size or
    mtime) is accepted only if every time indexed so far is unchanged (a
    running simulation appending dumps): before a cached index is used again,
    the address of its ``t/<time>`` header, a digest of the bytes there, the
    byte offsets of the first and last patch of every level and a fingerprint
    of their data are read again and compared. Any difference means the file
    was rewritten: PhareseeFileChanged, now and for every later call.
    """

    def __init__(self, path, *, nthreads=None, raw_io=None, prefetch="auto"):
        self._identity = file_identity(str(path))
        self._opened = self._identity  # the time addresses in memory describe this state
        self._path = self._identity[0]
        if nthreads is None:
            nthreads = rawio.io_threads()
        self._nthreads = max(1, int(nthreads))
        if raw_io is None:
            raw_io = rawio.raw_io_enabled()
        self._raw_io = bool(raw_io)
        if prefetch not in ("auto", "never", "walk"):
            raise ValueError(f"prefetch={prefetch!r}: use 'auto', 'never' or 'walk'")
        self._prefetch_mode = prefetch
        self._lock = threading.Lock()
        self._indexes = {}  # time key -> _CachedTime
        self._rewritten = None  # PhareseeFileChanged once a rewrite was detected
        self._prefetched = set()
        try:
            self._read_file_header()
        except (FormatError, PhareseeIOError):
            raise
        except (OSError, RuntimeError, KeyError) as e:
            # damaged metadata (e.g. a file truncated inside its headers) surfaces as
            # raw h5py errors while walking the groups
            raise PhareseeIOError(
                f"{self._path}: cannot read the file structure ({e}); the file is damaged "
                "or still being written"
            ) from e

    # ------------------------------------------------------------------ meta
    def _read_file_header(self):
        import h5py

        path = self._path
        with _open_h5(path) as f:
            attrs = {k: f.attrs[k] for k in f.attrs}
            for required in ("cell_width", "dimension", "domain_box"):
                if required not in attrs:
                    raise FormatError(
                        f"{path}: root attribute {required!r} missing; not a phareh5 file"
                    )
            layout = attrs.get("layoutType", "yee")
            if isinstance(layout, bytes):
                layout = layout.decode()
            if str(layout) != "yee":
                raise FormatError(f"{path}: layoutType={layout!r}, only 'yee' is supported")
            if "t" not in f:
                raise FormatError(f"{path}: no group /t; not a phareh5 field file")
            tg = h5py.h5g.open(f.id, b"t")
            names = []
            tg.links.iterate(names.append)
            keys = [n.decode() for n in names]
            addrs = [tg.links.get_info(n).u for n in names]
            try:
                values = np.array([float(k) for k in keys], dtype=np.float64)
            except ValueError as e:
                raise FormatError(f"{path}: time group name not a number ({e})") from None
            order = np.argsort(values, kind="stable")
            self._time_keys = [keys[i] for i in order]
            self._times = values[order]
            self._time_addr = {keys[i]: int(addrs[i]) for i in order}
            datasets, dtypes = self._first_patch_datasets(f)

        dim = int(np.asarray(attrs["dimension"]).ravel()[0])
        if dim not in (1, 2, 3):
            raise FormatError(f"{path}: dimension={dim}, expected 1, 2 or 3")
        dl0 = np.asarray(attrs["cell_width"], dtype=np.float64).reshape(-1)
        cells = np.asarray(attrs["domain_box"], dtype=np.int64).reshape(-1) + 1
        if dl0.size != dim or cells.size != dim:
            raise FormatError(f"{path}: cell_width/domain_box do not have {dim} values")
        periodic = np.asarray(attrs.get("periodicities", np.ones(dim)), dtype=np.int64)
        periodic = periodic.reshape(-1).astype(bool)
        interp = attrs.get("interpOrder")
        interp = None if interp is None else int(np.asarray(interp).ravel()[0])
        model = "hybrid" if interp is not None else "mhd"

        names, ordered = component_names(datasets)
        dtypes = tuple(dtypes[datasets.index(d)] for d in ordered)
        if len(set(dtypes)) != 1:
            raise FormatError(f"{path}: components have different dtypes {dtypes}")
        self._datasets = ordered
        self._dtypes = dtypes
        self._dual = tuple(dual_axes(model, d, dim) for d in ordered)

        # byte region of each time, bounded by consecutive t/<time> headers
        size = self._identity[2]
        sorted_addr = np.array(sorted(self._time_addr.values()) + [size], dtype=np.int64)
        self._region = {}
        for key, a in self._time_addr.items():
            nxt = sorted_addr[np.searchsorted(sorted_addr, a, side="right")]
            self._region[key] = (a, int(max(nxt, a)))

        known = {"cell_width", "dimension", "domain_box", "interpOrder", "layoutType", "origin",
                 "periodicities"}
        extra = {k: (v.decode() if isinstance(v, bytes) else v) for k, v in attrs.items()
                 if k not in known}
        quantity = os.path.splitext(os.path.basename(path))[0]
        dtype = np.dtype(dtypes[0]).newbyteorder("=")
        raw = "on" if self._raw_io else "off"
        signature = f"phareh5 {dtype.name} {len(ordered)} datasets/patch, raw io {raw}"
        self._meta = FileMeta(
            path=path,
            format="phareh5",
            quantity=quantity,
            dim=dim,
            dl0=dl0,
            domain_cells0=cells,
            periodic=periodic,
            interp_order=interp,
            model=model,
            names=names,
            dtype=dtype,
            colocated=True,
            writer_signature=signature,
            attrs=extra,
        )

    def _first_patch_datasets(self, f):
        """Dataset names and dtypes of the first patch that has datasets."""
        import h5py

        for key in self._time_keys:
            tg = h5py.h5g.open(f.id, f"t/{key}".encode())
            levels = []
            tg.links.iterate(levels.append)
            for ln in sorted(levels, key=_parse_level_name):
                lg = h5py.h5g.open(tg, ln)
                patches = []
                lg.links.iterate(patches.append)
                for pn in patches:
                    pg = h5py.h5g.open(lg, pn)
                    dnames = []
                    pg.links.iterate(dnames.append)
                    if dnames:
                        dsets = [n.decode() for n in dnames]
                        dtypes = [h5py.h5d.open(pg, n).dtype for n in dnames]
                        return dsets, dtypes
        raise FormatError(f"{self._path}: no patch with datasets in any time; empty file?")

    def meta(self) -> FileMeta:
        """FileMeta of the file (read at open, no I/O)."""
        return self._meta

    def times(self) -> np.ndarray:
        """Sorted float64 dump times (from the "%.10f" group names, at open)."""
        return self._times.copy()

    def spec(self) -> tuple:
        """Picklable recipe: ``open_backend(spec[1], **spec[2])`` reopens the file."""
        return ("phareh5", self._path, {"raw_io": self._raw_io, "prefetch": self._prefetch_mode})

    # ------------------------------------------------------------------ times
    def time_key(self, time) -> str:
        """Name of the time group for `time` (float, or the group name itself).

        Matches when |time - t| <= ``rawio.time_tolerance``: max(2 ulp of float64 at
        t, 5e-11), or 2 ulp of float32 when `time` is a float32 numpy value. Raises
        TimeNotFoundError listing the 3 nearest dump times otherwise.
        """
        if isinstance(time, (str, bytes)):
            text = time.decode() if isinstance(time, bytes) else time
            if text in self._time_addr:
                return text
            try:
                time = float(text)
            except ValueError:
                raise TimeNotFoundError(f"{self._path}: no time {text!r}") from None
        if self._times.size == 0:
            raise TimeNotFoundError(f"{self._path}: the file has no dump time")
        i = rawio.nearest_time(self._times, time, np.float64)
        if i is not None:
            return self._time_keys[i]
        t = float(time)
        near = self._times[np.argsort(np.abs(self._times - t), kind="stable")[:3]]
        raise TimeNotFoundError(
            f"{self._path}: no dump at time {t!r}; nearest times are "
            f"{', '.join(repr(float(v)) for v in np.sort(near))}"
        )

    # ------------------------------------------------------------------ index
    def index(self, time) -> tuple:
        """Tuple of LevelIndex (levels with patches, increasing level) at `time`.

        lower/upper are (npatch, dim) int64 inclusive cell boxes in FILE order
        (patch groups sorted by (rank, localid)); payload is an H5LevelTable.
        Uses, in order: the in-memory index, the on-disk cache (indexcache),
        a low-level walk of the time group (then stored in the cache).
        Raises TimeNotFoundError, FormatError (e.g. a time without level 0,
        written by per-level dumps: not supported by this reader).
        """
        key = self.time_key(time)
        levels = self._level_indexes(key)
        if 0 not in levels:
            raise FormatError(
                f"{self._path}: time {key} has no level 0 (a per-level dump); this reader "
                "needs level 0 at every time. Use Run(path, legacy=True) for this time."
            )
        return tuple(levels[l] for l in sorted(levels))

    def _current_identity(self):
        """Stat the file: PhareseeFileChanged if it was replaced, shrunk or found
        rewritten earlier; returns the current FileIdentity."""
        if self._rewritten is not None:
            raise self._rewritten
        try:
            return check_identity(self._identity)
        except PhareseeFileChanged as e:  # replaced or shrunk: final as well
            with self._lock:
                self._rewritten = e
                self._indexes.clear()
            raise

    def _mark_rewritten(self, what):
        """Forget every cached index and fail this and every later call."""
        err = rawio.rewritten(self._path, what)
        with self._lock:
            self._rewritten = err
            self._indexes.clear()
        raise err

    def _time_digest(self, key):
        """sha256 of the bytes at the t/<key> header, at most 4 KiB, never past the
        next time's header (as seen at open): bytes a later dump cannot touch."""
        start, stop = self._region[key]
        n = max(1, min(indexcache.HEADER_BYTES, stop - start))
        return indexcache.header_digest(self._path, start, n)

    def _check_time_address(self, key):
        """The t/<key> header is still where it was at open (else: rewritten)."""
        import h5py

        try:
            with _open_h5(self._path) as f:
                tg = h5py.h5g.open(f.id, b"t")
                addr = int(tg.links.get_info(key.encode()).u)
        except (KeyError, OSError, RuntimeError):
            addr = None
        if addr != self._time_addr[key]:
            self._mark_rewritten(f"the group t/{key} moved or disappeared")

    def _fingerprint(self, levels) -> bytes:
        """First bytes of the first and the last dataset (by file offset) of every level."""
        offsets = []
        for li in levels.values():
            off = li.payload.offset[li.payload.offset >= 0]
            if off.size:
                offsets += [int(off.min()), int(off.max())]
        if not offsets:
            return b""
        fd = rawio.open_readonly(self._path)
        try:
            return b"".join(os.pread(fd, _PRINT_BYTES, o) for o in offsets)
        finally:
            os.close(fd)

    # DESIGN-NOTE: "dataset offsets re-read" is done on a sample (first and last patch of
    # every level, all components) rather than by a full walk (seconds per time on 10^4
    # patches); the header digest and the data fingerprint cover the rest.
    def _same_offsets(self, key, levels) -> bool:
        """Byte offsets of the first and last patch of every level, read again from the file."""
        try:
            with _open_h5(self._path) as f:
                for lvl, li in levels.items():
                    t = li.payload
                    for i in {0, len(t.rank) - 1}:
                        for c, dname in enumerate(t.datasets):
                            ds = f[f"t/{key}/pl{lvl}/{t.patch_name(i)}/{dname}"]
                            off = ds.id.get_offset()
                            if (-1 if off is None else int(off)) != int(t.offset[c, i]):
                                return False
        except (KeyError, OSError, RuntimeError):
            return False
        return True

    def _level_indexes(self, key) -> dict:
        now = self._current_identity()
        sig = (now.st_size, now.st_mtime_ns)
        with self._lock:
            got = self._indexes.get(key)
        if got is not None and got.sig == sig:
            return got.levels
        modified = rawio.modified_since(self._opened, now)
        if modified:
            self._check_time_address(key)
        if got is not None:
            # the file was modified since this index was made: an append keeps it valid
            if (self._time_digest(key) != got.digest
                    or not self._same_offsets(key, got.levels)
                    or self._fingerprint(got.levels) != got.prints):
                self._mark_rewritten(f"the patches or the data of time {key} differ")
            entry = _CachedTime(got.levels, sig, got.digest, got.prints)
        else:
            levels = self._build_level_indexes(key)
            entry = _CachedTime(levels, sig, self._time_digest(key), self._fingerprint(levels))
        with self._lock:
            if self._rewritten is not None:
                raise self._rewritten
            self._indexes[key] = entry
            self._identity = now  # grown or modified, but what was used is unchanged
        return entry.levels

    def _build_level_indexes(self, key) -> dict:
        """{level: LevelIndex} of time `key`: from the on-disk cache or a walk."""
        addr = self._time_addr[key]
        inode = self._identity[1]
        # DESIGN-NOTE: the cache key also holds the address of the next time's
        # header (-1 for the last time). A file rewritten in place (same inode)
        # with other boxes can keep the same first 4 KiB at the time header
        # (old-style groups keep link names elsewhere; reproduced in
        # test_io_indexcache), but the size of its time regions changes.
        nxt = self._region[key][1] if self._region[key][1] < self._opened[2] else -1
        cache_key = f"{key}|next={nxt}"
        digest = None
        tables = None
        if indexcache.enabled():
            try:
                digest = indexcache.header_digest(self._path, addr)
                arrays = indexcache.load(self._path, inode, cache_key, addr, digest)
                if arrays is not None:
                    tables = self._tables_from_arrays(key, arrays)
            except (OSError, KeyError, ValueError, FormatError):
                tables = None
        if tables is None:
            if self._prefetch_mode == "walk":
                self.prefetch(key)
            tables = self._walk(key)
            if digest is not None:
                indexcache.store(
                    self._path, inode, cache_key, addr, digest, self._tables_to_arrays(tables)
                )
        levels = {}
        for lvl, (lower, upper, table) in tables.items():
            if lower.shape[0] > 0:
                levels[lvl] = LevelIndex(level=lvl, lower=lower, upper=upper, payload=table)
        return levels

    def _walk(self, key) -> dict:
        """Low-level walk of t/<key>: {level: (lower, upper, H5LevelTable)} in FILE order."""
        import h5py

        h5g, h5a, h5d = h5py.h5g, h5py.h5a, h5py.h5d
        STATS["walks"] += 1
        dim = self._meta.dim
        dsets = [d.encode() for d in self._datasets]
        nc = len(dsets)
        expected_dtypes = [np.dtype(t) for t in self._dtypes]
        out = {}
        with _open_h5(self._path) as f:
            tg = h5g.open(f.id, f"t/{key}".encode())
            lnames = []
            tg.links.iterate(lnames.append)
            for ln in lnames:
                lvl = _parse_level_name(ln)
                lg = h5g.open(tg, ln)
                pnames = []
                lg.links.iterate(pnames.append)
                n = len(pnames)
                rank = np.empty(n, np.int64)
                local = np.empty(n, np.int64)
                lower = np.empty((n, dim), np.int64)
                upper = np.empty((n, dim), np.int64)
                ghosts = np.zeros((nc, n), np.int64)
                shape = np.empty((nc, n, dim), np.int64)
                offset = np.full((nc, n), -1, np.int64)
                gbuf = np.empty((), np.int64)
                for i, pn in enumerate(pnames):
                    rank[i], local[i] = _parse_patch_name(pn)
                    pg = h5g.open(lg, pn)
                    self._read_box_attr(pg, b"lower", lower[i], key, pn)
                    self._read_box_attr(pg, b"upper", upper[i], key, pn)
                    for c, dn in enumerate(dsets):
                        try:
                            d = h5d.open(pg, dn)
                        except KeyError:
                            # DESIGN-NOTE: master keeps patch groups without
                            # datasets as empty patches; a field file never has
                            # them, so here it is a FormatError
                            raise FormatError(
                                f"{self._path}: t/{key}/{ln.decode()}/{pn.decode()} has no dataset "
                                f"{dn.decode()} (other patches have {self._datasets})"
                            ) from None
                        dshape = d.shape
                        if len(dshape) != dim:
                            raise FormatError(
                                f"{self._path}: dataset {dn.decode()} of {pn.decode()} has "
                                f"shape {dshape}, expected {dim} dimensions"
                            )
                        shape[c, i] = dshape
                        if d.dtype != expected_dtypes[c]:
                            raise FormatError(
                                f"{self._path}: dataset {dn.decode()} of {pn.decode()} is "
                                f"{d.dtype}, other patches are {expected_dtypes[c]}"
                            )
                        o = d.get_offset()
                        if o is not None:
                            offset[c, i] = o
                        try:
                            a = h5a.open(d, b"ghosts")
                        except KeyError:
                            continue
                        a.read(gbuf)
                        ghosts[c, i] = int(gbuf)
                order = np.lexsort((local, rank))
                itemsize = np.array([t.itemsize for t in expected_dtypes], np.int64)
                shape = shape[:, order]
                table = H5LevelTable(
                    time_key=key,
                    level=lvl,
                    rank=rank[order],
                    localid=local[order],
                    ghosts=ghosts[:, order],
                    shape=shape,
                    offset=offset[:, order],
                    nbytes=np.prod(shape, axis=2) * itemsize[:, None],
                    datasets=tuple(self._datasets),
                    dtypes=tuple(expected_dtypes),
                )
                lo, up = lower[order], upper[order]
                self._check_shapes(table, lo, up)
                out[lvl] = (lo, up, table)
        return out

    def _read_box_attr(self, group, name, dest, key, pn):
        import h5py

        try:
            a = h5py.h5a.open(group, name)
        except KeyError:
            raise FormatError(
                f"{self._path}: patch t/{key}/.../{pn.decode()} has no attribute {name.decode()}"
            ) from None
        if a.dtype.kind in "iu" and a.shape == dest.shape:
            a.read(dest)
            return
        # older writers stored boxes as text such as "[0 0]"
        buf = np.empty(a.shape, dtype=a.dtype)
        a.read(buf)
        val = buf[()] if buf.shape == () else buf
        if isinstance(val, bytes):
            val = val.decode()
        if isinstance(val, str):
            val = [int(v) for v in val.strip("[]() ").replace(",", " ").split()]
        dest[...] = np.asarray(val, dtype=np.int64).reshape(dest.shape)

    def _check_shapes(self, table, lower, upper):
        """FormatError if a dataset's shape does not match its box, centering and ghosts."""
        nn = upper - lower + 2
        for c, dname in enumerate(table.datasets):
            dual = np.array(self._dual[c], dtype=np.int64)
            expect = nn - dual[None, :] + 2 * table.ghosts[c][:, None]
            bad = np.flatnonzero(np.any(table.shape[c] != expect, axis=1))
            if bad.size:
                i = int(bad[0])
                raise FormatError(
                    f"{self._path}: t/{table.time_key}/pl{table.level}/"
                    f"{table.patch_name(i)}/{dname}"
                    f" has shape {tuple(table.shape[c, i])}, expected {tuple(expect[i])} for box "
                    f"{tuple(lower[i])}..{tuple(upper[i])} and ghosts {int(table.ghosts[c, i])}"
                )

    def _tables_to_arrays(self, tables) -> dict:
        arrays = {
            "levels": np.array(sorted(tables), dtype=np.int64),
            "datasets": np.array(self._datasets),
            "dtypes": np.array([np.dtype(t).str for t in self._dtypes]),
        }
        for lvl, (lower, upper, t) in tables.items():
            p = f"l{lvl}_"
            arrays[p + "lower"] = lower
            arrays[p + "upper"] = upper
            arrays[p + "rank"] = t.rank
            arrays[p + "localid"] = t.localid
            arrays[p + "ghosts"] = t.ghosts
            arrays[p + "shape"] = t.shape
            arrays[p + "offset"] = t.offset
        return arrays

    def _tables_from_arrays(self, key, arrays) -> dict:
        if tuple(str(d) for d in arrays["datasets"]) != tuple(self._datasets):
            return None
        dtypes = tuple(np.dtype(str(t)) for t in arrays["dtypes"])
        if dtypes != tuple(np.dtype(t) for t in self._dtypes):
            return None
        itemsize = np.array([t.itemsize for t in dtypes], np.int64)
        out = {}
        for lvl in arrays["levels"].tolist():
            p = f"l{lvl}_"
            shape = arrays[p + "shape"].astype(np.int64)
            table = H5LevelTable(
                time_key=key,
                level=int(lvl),
                rank=arrays[p + "rank"].astype(np.int64),
                localid=arrays[p + "localid"].astype(np.int64),
                ghosts=arrays[p + "ghosts"].astype(np.int64),
                shape=shape,
                offset=arrays[p + "offset"].astype(np.int64),
                nbytes=np.prod(shape, axis=2) * itemsize[:, None],
                datasets=tuple(self._datasets),
                dtypes=dtypes,
            )
            lower = arrays[p + "lower"].astype(np.int64)
            upper = arrays[p + "upper"].astype(np.int64)
            self._check_shapes(table, lower, upper)
            out[int(lvl)] = (lower, upper, table)
        return out

    def _level(self, time, level):
        key = self.time_key(time)
        levels = self._level_indexes(key)
        if 0 not in levels:
            self.index(key)  # raises the FormatError
        if level not in levels:
            raise ValueError(
                f"{self._path}: no level {level} at time {key}; levels are {sorted(levels)}"
            )
        return key, levels[level]

    # ------------------------------------------------------------------ sizes
    def nbytes(self, time, level, file_ids, names=None) -> int:
        """Bytes read_into reads from the file for these patches (datasets with
        their ghosts) and components `names` (None: all). No I/O beyond the index."""
        # DESIGN-NOTE: bytes read from disk (Yee + ghosts), slightly more than
        # the node values written into `out`; staging stays <= 2 x BATCH_BYTES.
        _, li = self._level(time, level)
        ids = np.asarray(file_ids, dtype=np.int64)
        comps = self._component_rows(names)
        return int(li.payload.nbytes[comps][:, ids].sum())

    def _component_rows(self, names):
        """Sorted dataset rows of components `names` (None: all). Raises ValueError for
        an empty selection, KeyError for an unknown name."""
        all_names = self._meta.names
        if names is None:
            return list(range(len(all_names)))
        names = [names] if isinstance(names, str) else list(names)
        unknown = [n for n in names if n not in all_names]
        if unknown:
            raise KeyError(f"unknown components {unknown}; this file has {all_names}")
        if not names:
            raise ValueError("no component asked for")
        return sorted({all_names.index(n) for n in names})

    def reads_whole_level(self, time, level, file_ids) -> bool:
        """Always False: this reader reads exactly the datasets asked for (see
        ``VtkhdfBackend.reads_whole_level``)."""
        return False

    def time_nbytes(self, time) -> int:
        """Bytes of all datasets of all levels at `time` (from the index)."""
        key = self.time_key(time)
        return int(sum(li.payload.nbytes.sum() for li in self._level_indexes(key).values()))

    # ------------------------------------------------------------------ prefetch
    def prefetch(self, time) -> float:
        """Warm the page cache over the byte region of one time; returns the seconds spent.

        The region starts at the t/<time> object header and ends at the next
        time's header (or the end of the file). It is read with threaded
        preadv through one scratch window of at most 256 MiB, then discarded:
        results never depend on it.
        """
        import time as _time

        key = self.time_key(time)
        start, stop = self._region[key]
        if stop <= start:
            return 0.0
        t0 = _time.perf_counter()
        window = min(PREFETCH_WINDOW, stop - start)
        scratch = np.empty(window, dtype=np.uint8)
        fd = rawio.open_readonly(self._path)
        try:
            stop = min(stop, os.fstat(fd).st_size)
            with ThreadPoolExecutor(self._nthreads) as pool:
                for w0 in range(start, stop, window):
                    n = min(stop, w0 + window) - w0
                    rawio.read_ranges(fd, [w0], [n], scratch[:n], executor=pool)
        finally:
            os.close(fd)
        with self._lock:
            self._prefetched.add(key)
        return _time.perf_counter() - t0

    # ------------------------------------------------------------------ reads
    def _validate_request(self, li, file_ids, dst_start, out, names):
        ids = np.asarray(file_ids, dtype=np.int64).reshape(-1)
        dst = np.asarray(dst_start, dtype=np.int64).reshape(-1)
        npatch = li.lower.shape[0]
        if ids.size != dst.size:
            raise ValueError(f"file_ids ({ids.size}) and dst_start ({dst.size}) differ in length")
        if ids.size and (ids.min() < 0 or ids.max() >= npatch):
            raise IndexError(f"file_ids out of range [0, {npatch}) at level {li.level}")
        missing = [n for n in names if n not in out]
        if missing:
            raise KeyError(f"out has no buffer for components {missing}; expected {names}")
        nn = li.upper[ids] - li.lower[ids] + 2
        need = int((dst + np.prod(nn, axis=1)).max(initial=0))
        for n in names:
            buf = out[n]
            if buf.ndim != 1 or buf.size < need:
                raise ValueError(f"out[{n!r}] must be 1-D with at least {need} values")
        return ids, dst, nn

    def read_into(self, time, level, file_ids, out, dst_start) -> None:
        """Fill `out` with the node values of patches `file_ids` of `level` at `time`.

        file_ids  : int64 file-order patch indices (any order; sorted is fastest)
        out       : {component name: 1-D buffer} for a non-empty subset of
                    meta().names; only those components' datasets are read
        dst_start : int64, node offset in `out` of each patch
        Each patch is written x fastest (Fortran order), nodes lower..upper+1,
        no ghosts, in the dtype of `out` (Yee averages done in the file's
        dtype when equal). Reads the datasets at their byte offsets with
        threaded preadv (h5py read_direct where a dataset has no single offset).
        Raises PhareseeIOError (file, time, level, patch count) on read errors,
        PhareseeFileChanged if the file was replaced or truncated.
        """
        key, li = self._level(time, level)
        comps = self._component_rows(list(out))
        names = tuple(self._meta.names[c] for c in comps)
        ids, dst, nn = self._validate_request(li, file_ids, dst_start, out, names)
        table = li.payload
        for c in comps:
            dual = self._dual[c]
            if any(dual):
                g = table.ghosts[c, ids]
                if g.size and g.min() < 2:
                    raise FormatError(
                        f"{self._path}: dataset {table.datasets[c]} at t/{key}/pl{level} has "
                        f"{int(g.min())} ghost layers; at least 2 are needed to average a dual "
                        "direction to nodes without the outermost ghost layer"
                    )
        self._maybe_prefetch(key, table, ids)

        def convert(k, raw):
            n = nn[k]
            size = int(np.prod(n))
            for c, name in zip(comps, names):
                view = out[name][dst[k] : dst[k] + size].reshape(tuple(n), order="F")
                yee_to_nodes(raw[c], int(table.ghosts[c, ids[k]]), n, self._dual[c], view)

        self._read_patches(key, li, ids, comps, convert)

    def read_derived(self, kind, time, level, file_ids, out, dst_start) -> None:
        """Yee-exact derived quantity of an EM_B file, averaged to nodes, into `out`.

        kind "J_yee"    : master's GetJ stencil (pharesee/run/utils.py:6-115) on
                          the raw Yee B of each patch, then each J component
                          averaged to nodes; out = {"x","y","z"} (in 1D Jx = 0)
        kind "divB_yee" : master's divB stencil (run/utils.py:118-121) at cells,
                          averaged from the 2^dim cells around each node;
                          out = {"value"}. Equals the centred node divergence of
                          the node-averaged B (DESIGN section 2).
        Same arguments, layout and errors as read_into. Computed in float64
        (differences of B in the file's dtype, like master), written in out's
        dtype. Raises ValueError for another kind or a file that is not EM_B.
        """
        if kind not in ("J_yee", "divB_yee"):
            raise ValueError(f"kind={kind!r}: use 'J_yee' or 'divB_yee'")
        if tuple(self._datasets) != _B_DATASETS:
            raise ValueError(
                f"{self._path}: {kind} needs the magnetic field file (EM_B.h5), "
                f"this file has {self._datasets}"
            )
        key, li = self._level(time, level)
        names = ("x", "y", "z") if kind == "J_yee" else ("value",)
        ids, dst, nn = self._validate_request(li, file_ids, dst_start, out, names)
        table = li.payload
        dim = self._meta.dim
        dl = self._meta.dl0 / 2.0**level
        g_all = table.ghosts[:, ids]
        if g_all.size and (g_all.min() < 2 or np.any(g_all != g_all[0])):
            raise FormatError(
                f"{self._path}: {kind} needs the same ghost width >= 2 on Bx, By, Bz "
                f"(t/{key}/pl{level})"
            )
        self._maybe_prefetch(key, table, ids)
        j_dual = [tuple(d == a for d in range(dim)) for a in range(3)]
        all_dual = (True,) * dim

        def convert(k, raw):
            n = nn[k]
            size = int(np.prod(n))
            g = int(table.ghosts[0, ids[k]])
            views = [
                out[name][dst[k] : dst[k] + size].reshape(tuple(n), order="F") for name in names
            ]
            if kind == "J_yee":
                J = yee_current(raw, g, n - 1, dl)
                for a in range(3):
                    yee_to_nodes(J[a], g, n, j_dual[a], views[a])
            else:
                yee_to_nodes(yee_divergence(raw, dl), g, n, all_dual, views[0])

        comps = range(3) if kind == "J_yee" else range(dim)
        self._read_patches(key, li, ids, comps, convert)

    def _maybe_prefetch(self, key, table, ids):
        # DESIGN-NOTE: "whole-time requests (>= 50 % of the time's bytes)" is
        # judged per read_into call, i.e. per level (the store reads one level
        # at a time): it fires when one call asks for >= 50 % of all levels'
        # bytes. Callers that know they read a whole time call prefetch(time).
        if self._prefetch_mode == "never" or key in self._prefetched:
            return
        requested = table.nbytes[:, ids].sum()
        total = self.time_nbytes(key)
        if total and requested >= PREFETCH_FRACTION * total:
            self.prefetch(key)

    def _read_patches(self, key, li, ids, comps, convert):
        """Read the raw datasets `comps` of patches `ids` batch by batch and call
        convert(k, raw) for each patch k (raw: list indexed by component, raw
        arrays [x, y, z] with ghosts; None for components not read). Retries
        once when the file changes during the read."""
        for attempt in range(2):
            self._current_identity()
            before = os.stat(self._path)
            try:
                self._read_patches_once(key, li, ids, list(comps), convert)
            except (PhareseeFileChanged, FormatError):
                raise
            except OSError as e:
                raise PhareseeIOError(
                    f"{self._path}: read failed at time {key}, level {li.level}, "
                    f"{len(ids)} patches: {e}"
                ) from e
            after = os.stat(self._path)
            if after.st_ino == before.st_ino and after.st_size >= before.st_size:
                return
            if attempt == 1:
                check_identity(self._identity)  # raises when replaced or truncated
                raise PhareseeFileChanged(
                    f"{self._path} changed during two reads; create a new Run."
                )

    def _read_patches_once(self, key, li, ids, comps, convert):
        table = li.payload
        ncomp = len(table.datasets)
        if ids.size == 0:
            return
        nbytes = table.nbytes[comps][:, ids]  # (nc, n)
        offset = table.offset[comps][:, ids]
        per_patch = nbytes.sum(axis=0)
        bidx = (np.cumsum(per_patch) - 1) // BATCH_BYTES
        cuts = np.concatenate([[0], np.flatnonzero(np.diff(bidx)) + 1, [ids.size]])
        fd = rawio.open_readonly(self._path)
        try:
            # `coord` runs the reads of batch b+1 (fanned out on `pool`) while
            # this thread converts batch b; two staging buffers alternate
            with ThreadPoolExecutor(self._nthreads) as pool, ThreadPoolExecutor(1) as coord:
                buffers = [None, None]

                def submit(b):
                    return self._submit_batch(
                        coord, pool, fd, buffers, b % 2, offset, nbytes, cuts, b
                    )

                pending = submit(0)
                for b in range(len(cuts) - 1):
                    current = pending
                    if b + 2 < len(cuts):
                        pending = submit(b + 1)
                    buf, pos, future = current
                    future.result()
                    p0, p1 = cuts[b], cuts[b + 1]
                    self._read_fallbacks(key, table, ids, comps, offset, pos, buf, p0, p1)
                    for k in range(p0, p1):
                        raw = [None] * ncomp
                        i = ids[k]
                        for j, c in enumerate(comps):
                            raw[c] = np.ndarray(
                                tuple(table.shape[c, i]),
                                dtype=table.dtypes[c],
                                buffer=buf,
                                offset=int(pos[j, k - p0]),
                            )
                        convert(k, raw)
        finally:
            os.close(fd)

    def _submit_batch(self, coord, pool, fd, buffers, slot, offset, nbytes, cuts, b):
        """Start the reads of batch b into staging buffer `slot`; returns (buf, pos, future)."""
        p0, p1 = cuts[b], cuts[b + 1]
        off = offset[:, p0:p1]
        nb = nbytes[:, p0:p1]
        raw = (off >= 0) if self._raw_io else np.zeros(off.shape, bool)
        starts, stops, _, rpos, total = plan_coalesced(off[raw], nb[raw])
        fb_sizes = nb[~raw]
        size = total + int(fb_sizes.sum())
        buf = buffers[slot]
        if buf is None or buf.size < size:
            buf = np.empty(max(size, 1), dtype=np.uint8)
            buffers[slot] = buf
        pos = np.empty(off.shape, np.int64)
        pos[raw] = rpos
        pos[~raw] = total + np.cumsum(fb_sizes) - fb_sizes
        future = coord.submit(
            rawio.read_ranges, fd, starts, stops - starts, buf[:total],
            nthreads=self._nthreads, executor=pool,
        )
        return buf, pos, future

    def _read_fallbacks(self, key, table, ids, comps, offset, pos, buf, p0, p1):
        """h5py read_direct of the datasets of patches p0..p1 that were not read raw."""
        raw = (offset[:, p0:p1] >= 0) if self._raw_io else np.zeros((len(comps), p1 - p0), bool)
        todo = np.argwhere(~raw)
        if todo.size == 0:
            return
        with _open_h5(self._path) as f:
            for j, kk in todo:
                c = comps[j]
                i = ids[p0 + kk]
                path = f"t/{key}/pl{table.level}/{table.patch_name(i)}/{table.datasets[c]}"
                arr = np.ndarray(
                    tuple(table.shape[c, i]), dtype=table.dtypes[c], buffer=buf,
                    offset=int(pos[j, kk]),
                )
                f[path].read_direct(arr)
