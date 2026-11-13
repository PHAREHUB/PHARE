"""Flat per-level buffers: lazy allocation and on-demand filling from a file backend.

``LevelData`` holds one quantity (all its components) on one level at one time: one flat
array per component, laid out by a ``LevelGeometry`` (store order, x fastest per patch).
File-backed data is read patch by patch when first needed and then kept; arrays handed to
users are read-only. In-memory data (results of operations) is plain numpy and writeable.
"""

from __future__ import annotations

import contextlib
import mmap
import sys
import threading
from typing import NamedTuple

import numpy as np

from .plan import default_budget

#: ``LevelData.ensure`` reads every patch of the LevelData instead of the ones asked for
#: when the backend reads the whole level span for the request anyway
#: (``reads_whole_level``) and the whole LevelData takes at most this fraction of the
#: memory available: the second access to the level then reads nothing.
EXPAND_FRACTION = 0.25


def lazy_empty(n, dtype):
    """Uninitialised 1D array of ``n`` values whose memory is only used where written.

    Linux: an anonymous private mmap advised ``MADV_NOHUGEPAGE`` and wrapped by
    ``np.frombuffer`` (the mmap stays alive through the array's ``base``). Untouched pages
    cost no resident memory; a written patch costs about its bytes (with transparent huge
    pages, ``np.empty`` would make every scattered patch cost up to 2 MiB). The address
    space reserved still counts toward ``RLIMIT_AS`` / virtual size.
    Elsewhere: ``np.empty``. Raises MemoryError if the address space cannot be reserved.
    """
    dtype = np.dtype(dtype)
    n = int(n)
    if n == 0 or not sys.platform.startswith("linux"):
        return np.empty(n, dtype=dtype)
    nbytes = n * dtype.itemsize
    try:
        buf = mmap.mmap(-1, nbytes, flags=mmap.MAP_PRIVATE | mmap.MAP_ANONYMOUS)
    except OSError as e:
        raise MemoryError(f"cannot reserve {nbytes} bytes of address space: {e}") from e
    if hasattr(mmap, "MADV_NOHUGEPAGE"):
        with contextlib.suppress(OSError):  # kernel without THP: nothing to disable
            buf.madvise(mmap.MADV_NOHUGEPAGE)
    return np.frombuffer(buf, dtype=dtype, count=n)


def _readonly_view(a):
    v = a.view()
    v.flags.writeable = False
    return v


def _reopen(spec, time, level, names, file_ids=None, derived=None, dtype=None):
    """Rebuild a file-backed LevelData from ``LevelData.__reduce__`` output (no data)."""
    from .geometry import default_cache

    backend = open_from_spec(spec)
    meta = backend.meta()
    index = next((li for li in backend.index(time) if li.level == level), None)
    if index is None:
        raise ValueError(f"{meta.path}: no level {level} at time {time}")
    lower, upper = index.lower, index.upper
    if file_ids is not None:
        lower, upper = lower[file_ids], upper[file_ids]
    r = 2**level
    geom = default_cache().get(
        level, lower, upper, np.asarray(meta.dl0) / r,
        np.asarray(meta.domain_cells0) * r, meta.periodic,
    )
    return LevelData(
        geom, time, names, meta.dtype if dtype is None else dtype, backend=backend,
        file_ids=file_ids, derived=derived,
    )


# formats whose backends return ``(format, path, kwargs)`` from ``spec()``
SPEC_FORMATS = ("vtkhdf", "phareh5")


def open_from_spec(spec):
    """Reopen a backend from ``Backend.spec()``.

    Accepted forms:
      ``(format, path, kwargs)`` with format in ``SPEC_FORMATS`` (what ``VtkhdfBackend`` and
          ``H5Backend`` return) -> ``io.open_backend(path, **kwargs)``;
      ``(callable, *args)`` -> ``callable(*args)`` (test backends);
      ``(path, kwargs)`` / ``(path,)`` -> ``io.open_backend(path, **kwargs)``.
    Raises ValueError for anything else; the backend's own errors (file missing, file
    replaced since the spec was made) propagate.
    """
    # DESIGN-NOTE: DESIGN 5 only says spec() is "a picklable recipe to reopen"; these are
    # the forms the backends produce.
    spec = tuple(spec)
    if not spec:
        raise ValueError("empty backend spec: expected (format, path, kwargs) from spec()")
    if callable(spec[0]):
        return spec[0](*spec[1:])
    from .io import open_backend

    if spec[0] in SPEC_FORMATS and len(spec) == 3 and isinstance(spec[2], dict):
        return open_backend(spec[1], **spec[2])
    if len(spec) <= 2 and (len(spec) == 1 or isinstance(spec[1], dict)):
        return open_backend(spec[0], **(spec[1] if len(spec) == 2 else {}))
    raise ValueError(
        f"backend spec {spec!r} not understood: expected (format, path, kwargs) with format "
        f"in {SPEC_FORMATS}, (callable, *args) or (path, kwargs)"
    )


class LevelDataBase:
    """Reading, buffers and views of one LevelData: everything but its construction.

    Works on the state ``LevelData.__init__`` sets (geom, time, names, dtype, backend, budget,
    file_ids, derived, _lock, _have, _buf) and has no ``__init__``: ``LevelData`` owns that
    state, ``field.ComponentLevel`` (one component of a LevelData) delegates it to its parent.
    """

    @property
    def level(self):
        return self.geom.level

    @property
    def file_backed(self):
        return self.backend is not None

    # ---- state
    @property
    def _loaded(self):
        """bool ``(npatch,)`` (a new array): patches holding every component."""
        return self._have.all(axis=0)

    @property
    def loaded(self):
        """Read-only bool ``(npatch,)``, store order: which patches hold every component."""
        return _readonly_view(self._loaded)

    def loaded_nbytes(self):
        """Bytes of the loaded (component, patch) pairs."""
        size = self.geom.size
        return int(sum(int(size[row].sum()) for row in self._have)) * self.dtype.itemsize

    def nbytes(self):
        """Logical size of the whole level (all components), loaded or not."""
        return self.geom.nnodes * self.dtype.itemsize * len(self.names)

    # ---- reading
    def ensure(self, slots=None, names=None, slab=None, force=False):
        """Make sure the patches ``slots`` (store order; None = all) hold components ``names``.

        names : None (every component) or a component name / list of names: only those are
                read (a derived quantity always reads every component)
        slab  : reserved for partial reads of 3D patches; must be None
        force : go ahead even if the memory check says it does not fit
        Reads only what is missing: the patches lacking one of the wanted components, and
        for them the wanted components missing at one of them at least (patches whose
        component was already loaded get the same values again). One ``backend.read_into``
        call with sorted file-order ids and ``out`` holding exactly those components, after
        ``MemoryBudget.check``. Exception: when the backend reads the whole level span for
        this request anyway (``backend.reads_whole_level``) and those components on every
        patch fit in ``EXPAND_FRACTION`` (25 %) of the available memory, every patch of the
        LevelData missing one of them is read in that call (a later access to another
        region of the level then reads nothing). The loaded state is set only once the
        data is in place: on failure it is unchanged and the error propagates
        (PhareseeIOError from the backend, PhareseeMemoryError from the budget).
        """
        if slab is not None:
            raise NotImplementedError(
                "ensure(slab=...) is reserved for later; read whole patches (slab=None)"
            )
        rows = self._rows(names)
        if self.backend is None:
            return
        slots = self._slots(slots)
        if self._have[np.ix_(rows, slots)].all():
            return
        with self._lock:
            missing, rows = self._missing_rows(slots, rows)
            if len(missing) == 0:  # another thread read them meanwhile
                return
            budget = default_budget() if self.budget is None else self.budget
            file_ids, missing = self._file_order(missing)
            if self._expand_to_whole(file_ids, rows, budget):
                every = np.arange(self.geom.npatch, dtype=np.int64)
                file_ids, missing = self._file_order(self._missing_rows(every, rows)[0])
            needed = self._backend_nbytes(file_ids, rows)
            budget.check(
                needed,
                f"reading {len(file_ids)} patches of level {self.geom.level} at t={self.time:g}",
                force=force,
            )
            if self._buf is None:
                # only address space is reserved here; pages are used where patches land
                self._buf = {n: lazy_empty(self.geom.nnodes, self.dtype) for n in self.names}
            dst_start = np.ascontiguousarray(self.geom.offset[missing])
            if self.derived is None:
                out = {self.names[r]: self._buf[self.names[r]] for r in rows}
                self.backend.read_into(self.time, self.geom.level, file_ids, out, dst_start)
            else:
                self.backend.read_derived(
                    self.derived, self.time, self.geom.level, file_ids, self._buf, dst_start
                )
            self._have[np.ix_(rows, missing)] = True

    def missing(self, slots=None, names=None):
        """(slots, names) that ``ensure(slots, names)`` would read, before any widening
        to the whole level: sorted store slots and component names. No I/O."""
        rows = self._rows(names)
        if self.backend is None:
            return np.empty(0, dtype=np.int64), ()
        slots, rows = self._missing_rows(self._slots(slots), rows)
        return slots, tuple(self.names[r] for r in rows)

    def read_nbytes(self, slots=None, names=None):
        """Bytes ``ensure(slots, names)`` reads (without widening to the whole level); 0 for
        in-memory data or when everything is loaded. No I/O."""
        slots, names = self.missing(slots, names)
        if len(slots) == 0:
            return 0
        return self._backend_nbytes(self._file_order(slots)[0], self._rows(names))

    def _rows(self, names):
        """Sorted component rows wanted by ``names`` (all for a derived quantity). Raises
        KeyError for an unknown name."""
        if names is not None:
            wanted = [names] if isinstance(names, str) else list(names)
            unknown = [n for n in wanted if n not in self.names]
            if unknown:
                raise KeyError(f"unknown component {unknown}; this data has {self.names}")
        if names is None or self.derived is not None:
            return list(range(len(self.names)))
        return sorted({self.names.index(n) for n in wanted})

    def _slots(self, slots):
        if slots is None:
            return np.arange(self.geom.npatch, dtype=np.int64)
        slots = np.unique(np.asarray(slots, dtype=np.int64).ravel())
        if len(slots) and (slots[0] < 0 or slots[-1] >= self.geom.npatch):
            raise IndexError(f"patch slots must be in [0, {self.geom.npatch})")
        return slots

    def _missing_rows(self, slots, rows):
        """(slots lacking one of `rows`, the rows of `rows` lacking at one of them)."""
        have = self._have[np.ix_(rows, slots)]
        lack = ~have.all(axis=0)
        missing = slots[lack]
        if len(missing) == 0:
            return missing, []
        return missing, [r for r, h in zip(rows, have[:, lack]) if not h.all()]

    def _backend_nbytes(self, file_ids, rows):
        """``backend.nbytes`` of components ``rows`` (no ``names`` argument when all)."""
        if self.derived is not None or len(rows) == len(self.names):
            return self.backend.nbytes(self.time, self.geom.level, file_ids)
        names = tuple(self.names[r] for r in rows)
        return self.backend.nbytes(self.time, self.geom.level, file_ids, names=names)

    def _backend_ids(self, slots):
        """Backend file-order patch ids of store slots `slots`."""
        file_ids = self.geom.file_index[slots]
        return file_ids if self.file_ids is None else self.file_ids[file_ids]

    def _file_order(self, slots):
        """(backend file ids, slots), both sorted by file id."""
        file_ids = self._backend_ids(slots)
        order = np.argsort(file_ids)
        return file_ids[order], slots[order]

    def _expand_to_whole(self, file_ids, rows, budget):
        """True when the request should be widened to every patch of this LevelData.

        That is when the backend says it reads the whole level span for this request
        anyway (``backend.reads_whole_level``, e.g. the vtkhdf S1 strategy: the bytes
        between the first and the last requested patch are read, most of them thrown
        away) AND the components ``rows`` on all the patches of this LevelData take at most
        ``EXPAND_FRACTION`` of the memory available (``budget.available()``). Backends
        without the method, and derived quantities, are never expanded.
        """
        if self.derived is not None or len(file_ids) >= self.geom.npatch:
            return False
        whole = getattr(self.backend, "reads_whole_level", None)
        available = getattr(budget, "available", None)
        if whole is None or available is None:
            return False
        if not whole(self.time, self.geom.level, file_ids):
            return False
        all_ids = np.sort(self._backend_ids(np.arange(self.geom.npatch, dtype=np.int64)))
        total = self._backend_nbytes(all_ids, rows)
        return total <= EXPAND_FRACTION * available()

    def raw(self, name):
        """Flat array ``(nnodes,)`` of component ``name`` as it is now, WITHOUT reading.

        Its values are only valid at the patches where ``loaded`` is True (call ``ensure``
        on the slots you need first); elsewhere they are undefined. Read-only when
        file-backed, the array itself (writeable) when in memory. If nothing was read yet,
        the lazy buffers are allocated here (address space only, see ``lazy_empty``) so that
        later ``ensure`` calls fill this same array. Raises KeyError for an unknown name.
        """
        self._check_name(name)
        if self._buf is None:
            with self._lock:
                if self._buf is None:
                    self._buf = {n: lazy_empty(self.geom.nnodes, self.dtype) for n in self.names}
        b = self._buf[name]
        return _readonly_view(b) if self.backend is not None else b

    def buffer(self, name):
        """Whole flat array ``(nnodes,)`` of component ``name`` (reads the missing patches of
        that component first).

        Read-only when file-backed; the array itself (writeable) when in memory.
        """
        self._check_name(name)
        self.ensure(None, name)
        b = self._buf[name]
        return _readonly_view(b) if self.backend is not None else b

    def buffers(self, names=None):
        """``{name: buffer(name)}`` for ``names`` (None: every component, in ``names``
        order), the missing ones read together in one ``ensure`` call first. Use it rather
        than ``buffer`` in a loop when several components are needed: one read per level
        instead of one per component."""
        wanted = self.names if names is None else (
            (names,) if isinstance(names, str) else tuple(names))
        self.ensure(None, wanted)
        return {n: self.buffer(n) for n in wanted}

    def view(self, slot, name):
        """Patch ``slot`` (store order) of component ``name``, indexed ``[x, y, z]``, no copy.

        Reads that component of the patch first if needed; read-only when file-backed.
        """
        self._check_name(name)
        self.ensure([slot], name)
        b = self._buf[name]
        if self.backend is not None:
            b = _readonly_view(b)
        return self.geom.view(b, slot)

    # DESIGN-NOTE: unload on in-memory data does nothing (it could not be read back).
    def unload(self):
        """Forget the data of a file-backed LevelData (next access reads again).

        Drops the references to the buffers (views a user still holds stay valid). Does
        nothing on in-memory data, which could not be read back.
        """
        if self.backend is None:
            return
        with self._lock:
            self._buf = None
            self._have[:] = False

    def _check_name(self, name):
        if name not in self.names:
            raise KeyError(f"unknown component {name!r}; this data has {self.names}")

    # ---- pickling
    def __reduce__(self):
        if self.backend is not None:
            return (_reopen, (self.backend.spec(), self.time, self.geom.level, self.names,
                              self.file_ids, self.derived, self.dtype))
        return (LevelData.from_arrays, (self.geom, self.time, dict(self._buf)))

    def __repr__(self):
        kind = "file" if self.backend is not None else "memory"
        return (
            f"LevelData(level={self.geom.level}, t={self.time:g}, names={self.names}, "
            f"dtype={self.dtype}, {kind}, loaded {int(self._loaded.sum())}/{self.geom.npatch})"
        )


class LevelData(LevelDataBase):
    """One quantity (all components), one level, one time, as flat arrays in store order.

    geom    : LevelGeometry (layout of the buffers)
    time    : float, the dump time
    names   : tuple of component names (e.g. ("x", "y", "z"))
    dtype   : numpy dtype of the values
    backend : the file backend (None for in-memory data)
    budget  : MemoryBudget consulted before every read (default ``plan.default_budget()``)
    file_ids: None, or int64 ``(npatch,)``: the backend's file-order patch id of each patch
              of ``geom`` taken in ITS file order (``geom`` then describes a subset of the
              level's patches, e.g. those touching a box)
    derived : None, or the name of a derived quantity the backend computes while reading
              (``backend.read_derived(derived, ...)`` is called instead of ``read_into``)

    File-backed: buffers are reserved with ``lazy_empty`` on the first read; ``ensure`` reads
    the (component, patch) pairs not yet loaded (one ``backend.read_into`` call), ``loaded``
    says which patches hold every component; arrays returned by ``buffer`` and ``view`` are
    read-only. In-memory (``from_arrays``): everything loaded, arrays writeable.
    Thread-safe: concurrent ``ensure`` calls read every patch once.
    """

    def __init__(self, geom, time, names, dtype, backend=None, budget=None, buffers=None,
                 file_ids=None, derived=None):
        if file_ids is not None:
            file_ids = np.asarray(file_ids, dtype=np.int64)
            if file_ids.shape != (geom.npatch,):
                raise ValueError(f"file_ids must have one entry per patch ({geom.npatch})")
        self.file_ids = file_ids
        self.derived = derived
        self.geom = geom
        self.time = float(time)
        self.names = tuple(names)
        self.dtype = np.dtype(dtype)
        self.backend = backend
        self.budget = budget
        self._lock = threading.Lock()
        # one row per component (order of ``names``), one column per patch (store order)
        self._have = np.zeros((len(self.names), geom.npatch), dtype=bool)
        self._buf = None
        if buffers is not None:
            self._buf = dict(buffers)
            self._have[:] = True

    # ---- construction
    @classmethod
    def from_arrays(cls, geom, time, arrays):
        """In-memory LevelData from ``{name: flat array (geom.nnodes,)}`` (store order).

        The arrays are used as given (no copy) and stay writeable. Raises ValueError on a
        length mismatch or mixed dtypes.
        """
        arrays = {k: np.asarray(v) for k, v in arrays.items()}
        if not arrays:
            raise ValueError("from_arrays needs at least one component")
        for name, a in arrays.items():
            if a.ndim != 1 or len(a) != geom.nnodes:
                raise ValueError(
                    f"component {name!r}: expected a flat array of {geom.nnodes} nodes, "
                    f"got shape {a.shape}"
                )
        dtypes = {a.dtype for a in arrays.values()}
        if len(dtypes) != 1:
            raise ValueError(f"components have different dtypes {dtypes}")
        return cls(geom, time, tuple(arrays), dtypes.pop(), buffers=arrays)

    @classmethod
    def from_backend(cls, backend, geom, time, names=None, dtype=None, budget=None,
                     file_ids=None, derived=None):
        """File-backed LevelData for level ``geom.level`` of ``backend`` at ``time``. No I/O.

        ``file_ids`` and ``derived``: see the class docstring.
        """
        if names is None or dtype is None:
            meta = backend.meta()
            names = meta.names if names is None else names
            dtype = meta.dtype if dtype is None else dtype
        return cls(geom, time, names, dtype, backend=backend, budget=budget,
                   file_ids=file_ids, derived=derived)


class Snapshot(NamedTuple):
    """One quantity at one time: the unit that operations and views exchange.

    hier   : HierGeometry (levels present at that time + extent)
    levels : tuple of LevelData, one per entry of ``hier.levels``, in the same order
             (``levels[k].geom is hier.levels[k]``), all with the same component names.
    """

    hier: object
    levels: tuple

    @property
    def names(self):
        """Component names, e.g. ("value",) or ("x", "y", "z")."""
        return self.levels[0].names

    @property
    def time(self):
        """Time of the snapshot (float)."""
        return self.levels[0].time
