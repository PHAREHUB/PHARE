"""User-facing AMR fields: ``ScalarField``, ``VectorField``, ``TensorField`` (base ``AMRField``).

A field is one quantity (all its components) at one or several times. Each time is a
``store.Snapshot``: a ``HierGeometry`` (patch boxes of every level, extent) and one
``LevelData`` per level (flat buffers, store order, x fastest inside a patch). Fields built
from a file are lazy: construction reads only the small box tables (``backend.index``);
field bytes are read at first use, for the patches needed, then kept. Results of operations
are in-memory fields; they chain (numpy protocol in ``_protocol.py``).

Indexing convention everywhere: patch arrays are indexed ``[x, y, z]``, node ``i`` of level
``l`` sits at ``i * dl0 / 2**l`` (physical units of the simulation).
"""

from __future__ import annotations

import numbers

import numpy as np

from . import plan
from ._protocol import NumpyOps
from .boxes import BinIndex, PhysBox, box_flat_indices, periodic_shifts
from .compat import MasterAPI
from .errors import TimeNotFoundError
from .geometry import HierGeometry, default_cache
from .store import LevelData, Snapshot, open_from_spec

SCALAR_NAMES = ("value",)
VECTOR_NAMES = ("x", "y", "z")
TENSOR6_NAMES = ("xx", "xy", "xz", "yy", "yz", "zz")
TENSOR9_NAMES = ("xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy", "zz")

_MIN_TIME_TOLERANCE = 5e-11


# --------------------------------------------------------------------------- time matching


def time_tolerance(t, time_dtype=np.float64):
    """Matching tolerance at time ``t``: max(2 ulp of ``time_dtype`` at t, 5e-11).

    Same rule as the backends (``io.vtkhdf.time_tolerance``, ``io.phareh5``).
    """
    ulp = float(np.spacing(np.abs(np.asarray(t, dtype=np.dtype(time_dtype)))))
    return max(2.0 * ulp, _MIN_TIME_TOLERANCE)


def match_time(times, t, time_dtype=np.float64, where=""):
    """Index of the entry of ``times`` (float64 array) equal to ``t`` within ``time_tolerance``.

    Raises TimeNotFoundError (a KeyError and a ValueError) listing the 3 nearest times.
    """
    from .io.rawio import time_tolerance as reader_tolerance

    times = np.asarray(times, dtype=np.float64)
    given = t  # a numpy float32 time is matched at float32 precision, as the readers do
    try:
        t = float(t)
    except (TypeError, ValueError):
        raise TypeError(f"a time must be a number, got {t!r}") from None
    if len(times):
        i = int(np.argmin(np.abs(times - t)))
        if abs(times[i] - t) <= reader_tolerance(t, time_dtype, given=given):
            return i
    near = np.sort(times[np.argsort(np.abs(times - t), kind="stable")[:3]])
    near_txt = ", ".join(repr(float(v)) for v in near) or "none"
    prefix = f"{where}: " if where else ""
    raise TimeNotFoundError(
        f"{prefix}no data at time {t!r}; nearest times: {near_txt}. List them with times()."
    )


def _times_list(times):
    """``times`` argument -> list of floats (a number, or a sequence of numbers)."""
    if isinstance(times, (numbers.Number, np.generic)) or (
        isinstance(times, np.ndarray) and times.ndim == 0
    ):
        return [float(times)]
    if isinstance(times, (str, bytes)):
        raise TypeError(f"times={times!r}: give a number or a list of numbers")
    return [float(t) for t in times]


# --------------------------------------------------------------------------- spec reopening


def open_spec(spec):
    """Reopen a backend from ``backend.spec()``.

    Accepts ``(format, path, kwargs)`` (what the vtkhdf and phareh5 backends return) and the
    forms of ``store.open_from_spec`` (``(callable, *args)``, ``(path, kwargs)``).
    """
    # DESIGN-NOTE: store.open_from_spec does not understand the (format, path, kwargs) spec of
    # the real backends (it calls open_backend(format, **path)); handled here.
    if (len(spec) == 3 and isinstance(spec[0], str) and isinstance(spec[2], dict)
            and not callable(spec[0])):
        from .io import open_backend

        return open_backend(spec[1], **spec[2])
    return open_from_spec(spec)


# --------------------------------------------------------------------------- component views


class ComponentLevel(LevelData):
    """One component of a LevelData seen as a scalar LevelData ``("value",)``, sharing buffers.

    Nothing is copied and nothing is read when it is made. Reading goes through the parent
    and reads this component only (a derived parent reads all its components);
    ``loaded`` says which patches hold this component; ``buffer``/``view`` return the
    parent's arrays for that component (read-only when file-backed).
    """

    # DESIGN-NOTE: LevelData has no "one component, shared buffers, no read" form; this
    # subclass gives it. unload() unloads the parent (the buffers are shared).
    def __init__(self, parent, component):  # state lives in the parent: no super().__init__
        if component not in parent.names:
            raise KeyError(f"unknown component {component!r}; this data has {parent.names}")
        self.parent = parent
        self.component = component
        self.geom = parent.geom
        self.time = parent.time
        self.names = SCALAR_NAMES
        self.dtype = parent.dtype

    # state delegated to the parent
    backend = property(lambda self: self.parent.backend)
    budget = property(lambda self: self.parent.budget)
    file_ids = property(lambda self: self.parent.file_ids)
    derived = property(lambda self: self.parent.derived)
    _lock = property(lambda self: self.parent._lock)

    @property
    def _have(self):
        """The parent's loaded row of this component, as a (1, npatch) view."""
        i = self.parent.names.index(self.component)
        return self.parent._have[i : i + 1]

    @property
    def _buf(self):
        buf = self.parent._buf
        return None if buf is None else {"value": buf[self.component]}

    def nbytes(self):
        return self.geom.nnodes * self.dtype.itemsize

    def ensure(self, slots=None, names=None, slab=None, force=False):
        if names is not None:
            self._check_names(names)
        self.parent.ensure(slots, self.component, slab, force)

    def missing(self, slots=None, names=None):
        if names is not None:
            self._check_names(names)
        slots, comps = self.parent.missing(slots, self.component)
        return slots, (SCALAR_NAMES if comps else ())

    def read_nbytes(self, slots=None, names=None):
        if names is not None:
            self._check_names(names)
        return self.parent.read_nbytes(slots, self.component)

    def _check_names(self, names):
        bad = [n for n in ([names] if isinstance(names, str) else names) if n != "value"]
        if bad:
            raise KeyError(f"unknown component {bad}; this data has {SCALAR_NAMES}")

    def buffer(self, name):
        self._check_names(name)
        return self.parent.buffer(self.component)

    def view(self, slot, name):
        self._check_names(name)
        return self.parent.view(slot, self.component)

    def raw(self, name):
        self._check_names(name)
        return self.parent.raw(self.component)

    def unload(self):
        self.parent.unload()

    def __reduce__(self):
        return (ComponentLevel, (self.parent, self.component))

    def __repr__(self):
        return f"ComponentLevel({self.component!r} of {self.parent!r})"


def component_level(ld, name):
    """Scalar LevelData of component ``name`` of ``ld``, sharing its buffers, without reading."""
    if ld.backend is None and not isinstance(ld, ComponentLevel):
        ld._check_name(name)
        return LevelData.from_arrays(ld.geom, ld.time, {"value": ld.buffer(name)})
    return ComponentLevel(ld, name)


# --------------------------------------------------------------------------- master's selection


class MasterSelection:
    """Master's ``selection_box`` patch rule (fromh5.py:235-247), applied to box tables.

    lower, upper : inclusive integer bounds of master's ``Box`` (physical units; master's
                   ``Box`` truncates float bounds to int, core/box.py:15-16, issue #1451).
    A patch of a level of mesh size ``dl`` is selected when the integer box
    ``Box(origin, origin + nbrCells * dl)`` (``origin = lower * dl``, truncated to int the same
    way) intersects ``[lower, upper]`` inclusively, as master does. So float bounds are
    truncated and patches that only touch the box are kept. Picklable, no I/O.
    """

    __slots__ = ("lower", "upper")

    def __init__(self, lower, upper):
        self.lower = tuple(int(v) for v in lower)
        self.upper = tuple(int(v) for v in upper)

    def selects(self, lower, upper, dl):
        """Bool ``(npatch,)``: which patches (int64 ``(npatch, dim)`` cell boxes) master keeps."""
        lower = np.asarray(lower, dtype=np.int64)
        upper = np.asarray(upper, dtype=np.int64)
        dl = np.asarray(dl, dtype=np.float64)
        origin = lower * dl
        pos_upper = origin + (upper - lower + 1) * dl
        # master: Box(origin, pos_upper) casts with astype(int) (truncation toward zero)
        plo, phi = origin.astype(int), pos_upper.astype(int)
        lo = np.maximum(plo, np.asarray(self.lower))
        hi = np.minimum(phi, np.asarray(self.upper))
        return (lo <= hi).all(axis=1)

    def __getstate__(self):
        return (self.lower, self.upper)

    def __setstate__(self, state):
        self.lower, self.upper = state

    def __eq__(self, other):
        return (isinstance(other, MasterSelection) and self.lower == other.lower
                and self.upper == other.upper)

    def __hash__(self):
        return hash((self.lower, self.upper))

    def __repr__(self):
        return f"MasterSelection(lower={self.lower}, upper={self.upper})"


# --------------------------------------------------------------------------- pickling helpers


def _rebuild_from_spec(spec, kwargs, component):
    field = AMRField.from_backend(open_spec(spec), **kwargs)
    return field[component] if component is not None else field


def _rebuild_from_levels(cls_name, snap_states, attrs):
    snaps = []
    for extent, levels in snap_states:
        hier = HierGeometry(tuple(ld.geom for ld in levels), extent)
        snaps.append(Snapshot(hier, tuple(levels)))
    field = make_field(snaps, **attrs)
    if type(field).__name__ != cls_name:  # pragma: no cover - make_field is deterministic
        raise TypeError(f"unpickled a {type(field).__name__}, expected {cls_name}")
    return field


# --------------------------------------------------------------------------- the field


def _fmt_int(n):
    return f"{int(n):,}".replace(",", " ")


def _fmt_bytes(n):
    return plan.format_bytes(int(n))


def _fmt_vec(v):
    return "(" + ", ".join(f"{float(x):g}" for x in np.asarray(v).ravel()) + ")"


class AMRField(NumpyOps, MasterAPI):
    """One quantity (all components) at one or several times, on an AMR hierarchy.

    Use ``ScalarField``, ``VectorField``, ``TensorField`` (``make_field`` picks the class
    from the component names); build from a file with ``AMRField.from_backend``.

    Attributes
      names        component names: ("value",), ("x", "y", "z"), tensor ("xx", "xy", ...)
      quantity     quantity name (e.g. "EM_B"; "EM_B_x" for a component) or None
      meta         io.FileMeta of the source file, or None
      path, format source file path and format ("vtkhdf", "phareh5"), or None
      stencil      "node" / "yee" for derived quantities, else None
      extent       PhysBox restricting the field, or None (whole domain)
      ndim, dl     dimension; level-0 mesh size (dim,) float64
      domain       physical domain size (dim,) float64; periodic (dim,) bool
      interp_order int, or None when the file does not say (MHD)
      sim          simulation object, or None (always None on vtkhdf)
    Data is node centred, no ghost nodes, ``[x, y, z]`` indexing; file-backed arrays are
    read-only, ``copy()`` gives a writeable in-memory field.
    Equality ``==`` between fields is identity (as on master) and fields hash by identity;
    ``field == number`` / ``field == array`` raises TypeError: elementwise comparison is
    ``np.equal`` / ``np.not_equal`` / ``np.isclose`` (fields of bool).
    Arithmetic between two times (``B1 - B0``) is done on the common refinement of the two
    layouts (cell-box intersections of their patches, DESIGN 2), so ``(B1 - B0).integral()``
    integrates over that common layout and can differ from ``B1.integral() - B0.integral()``
    where the refined regions moved (each of those uses its own best data).
    """

    _expected = None  # names check of subclasses: callable(names) -> bool, and a description

    def __init__(self, snaps, *, meta=None, quantity=None, stencil=None, sim=None):
        snaps = tuple(snaps)
        if not snaps:
            raise ValueError(f"{type(self).__name__} needs at least one time (Snapshot)")
        for s in snaps:
            if not isinstance(s, Snapshot):
                raise TypeError(
                    f"{type(self).__name__} is built from store.Snapshot objects, got "
                    f"{type(s).__name__}; read files with AMRField.from_backend"
                )
        names = tuple(snaps[0].names)
        if any(tuple(s.names) != names for s in snaps[1:]):
            raise ValueError("every time of a field must have the same component names")
        self._check_names(names)
        times = np.array([s.time for s in snaps], dtype=np.float64)
        if len(set(times.tolist())) != len(times):
            raise ValueError(f"a time appears twice in {times.tolist()}; give each time once")
        self._snaps = snaps
        self._times = times
        self._times.flags.writeable = False
        self.meta = meta
        self.quantity = quantity
        self.stencil = stencil
        self.sim = sim
        self._time_dtype = np.dtype(np.float64)
        self._origin = None  # how to reopen this field from its file (pickling), or None
        # name of this component in its parent ("x" for B.x), or None: master's views of a
        # component answer the parent's key (compat.py); results of operations have none
        self._component_alias = None

    @classmethod
    def _check_names(cls, names):
        if cls._expected is not None and not cls._expected[0](names):
            raise TypeError(
                f"a {cls.__name__} needs {cls._expected[1]}, got components {names}; use "
                "make_field(...) to get the class matching the components"
            )

    # ------------------------------------------------------------------ construction
    @classmethod
    def from_backend(cls, backend, times=None, extent=None, cache=None, quantity=None,
                     derived=None, names=None, dtype=None, stencil=None, sim=None,
                     selection=None):
        """Field of a file backend at some times, reading NO field bytes.

        backend  : a Backend (``io.open_backend(path)``)
        times    : None = every time of the file (sorted); a number or a list of numbers, kept
                   in the caller's order; matched with the backends' tolerance
                   (TimeNotFoundError listing the nearest times otherwise)
        extent   : None, ``dict(x=(x0, x1), ...)`` or PhysBox (physical units). Each level then
                   keeps only the patches whose cell box meets the box's cells (PhysBox rule);
                   the first level left without patches ends the hierarchy.
        cache    : GeometryCache (default: the module-level one)
        quantity : name shown in repr (default: the file's quantity)
        derived, names, dtype : passed to ``LevelData.from_backend`` (backend-computed
                   quantities such as Yee J on phareh5)
        stencil, sim : stored on the field
        selection: None, or a ``MasterSelection`` (master's ``selection_box`` rule) choosing
                   the patches instead of the ``extent`` rule; ``extent`` stays the field's
                   extent
        Reads per time only the box tables (``backend.index``). The returned class follows
        the component names when called on AMRField, else must match ``cls`` (TypeError).
        """
        meta = backend.meta()
        cache = default_cache() if cache is None else cache
        # DESIGN-NOTE: the Backend protocol does not expose the dtype of the file's times;
        # vtkhdf keeps it as ``_time_dtype`` (float32 times on current files), used for the
        # tolerance of ``at(t)``; float64 otherwise.
        time_dtype = np.dtype(getattr(backend, "_time_dtype", np.float64))
        file_times = np.asarray(backend.times(), dtype=np.float64)
        where = getattr(meta, "path", "")
        if times is None:
            wanted = file_times.tolist()
        else:
            wanted = [float(file_times[match_time(file_times, t, time_dtype, where)])
                      for t in _times_list(times)]
        box = PhysBox.parse(extent, meta.dim)
        snaps = [_snapshot_from_backend(backend, meta, t, box, cache, names, dtype, derived,
                                        selection)
                 for t in wanted]
        kw = dict(meta=meta, quantity=meta.quantity if quantity is None else quantity,
                  stencil=stencil, sim=sim)
        field = make_field(snaps, **kw) if cls is AMRField else cls(snaps, **kw)
        field._time_dtype = time_dtype
        field._origin = dict(
            spec=backend.spec(),
            kwargs=dict(times=wanted, extent=box, quantity=field.quantity, derived=derived,
                        names=names, dtype=dtype, stencil=stencil, selection=selection),
            component=None,
        )
        return field

    def _sharing(self, snaps):
        """Same field (class, meta, quantity, stencil, sim, origin, component name) on other
        Snapshots of the same times, e.g. Snapshots another field already holds. No I/O."""
        snaps = tuple(snaps)
        if [float(s.time) for s in snaps] != [float(t) for t in self._times]:
            raise ValueError("_sharing: the snapshots must be at the field's times")
        field = self._like(snaps)
        field._origin = self._origin
        field._component_alias = self._component_alias
        return field

    def _like(self, snaps, **overrides):
        """New field (class chosen from the component names) with this field's meta,
        quantity, stencil and sim unless overridden."""
        kw = dict(meta=self.meta, quantity=self.quantity, stencil=self.stencil, sim=self.sim)
        kw.update(overrides)
        field = make_field(snaps, **kw)
        field._time_dtype = self._time_dtype
        return field

    # ------------------------------------------------------------------ times
    def times(self):
        """Times held by the field, float64 array, in the caller's order (a copy)."""
        return self._times.copy()

    def _index_of(self, time):
        return match_time(self._times, time, self._time_dtype, where=self._label())

    def snapshot(self, time=None):
        """The Snapshot at ``time``; with ``time=None`` the only one (ValueError when the field
        holds several times: the message lists them). TimeNotFoundError if absent."""
        if time is None:
            if len(self._snaps) == 1:
                return self._snaps[0]
            raise ValueError(
                f"{self._label()} holds {len(self._snaps)} times "
                f"{[float(t) for t in self._times]}; pass time=... or use field.at(t)"
            )
        return self._snaps[self._index_of(time)]

    def snapshots(self):
        """Tuple of Snapshot, one per time, in the field's order."""
        return self._snaps

    def at(self, t):
        """The same field at one time ``t`` only (shares data, reads nothing)."""
        i = self._index_of(t)
        field = self._like([self._snaps[i]])
        field._component_alias = self._component_alias
        if self._origin is not None:
            field._origin = dict(self._origin, kwargs=dict(self._origin["kwargs"],
                                                           times=[float(self._times[i])]))
        return field

    def hier(self, time=None):
        """HierGeometry at ``time`` (single-time rule of ``snapshot``)."""
        return self.snapshot(time).hier

    def levels_data(self, time=None):
        """Tuple of LevelData (one per level present) at ``time`` (rule of ``snapshot``)."""
        return self.snapshot(time).levels

    # ------------------------------------------------------------------ description
    @property
    def names(self):
        return tuple(self._snaps[0].names)

    @property
    def dtype(self):
        """numpy dtype of the values."""
        return self._snaps[0].levels[0].dtype

    @property
    def path(self):
        return None if self.meta is None else self.meta.path

    @property
    def format(self):
        return None if self.meta is None else self.meta.format

    @property
    def extent(self):
        return self._snaps[0].hier.extent

    @property
    def ndim(self):
        return self._snaps[0].hier.dim

    @property
    def dl(self):
        """Level-0 mesh size, float64 (dim,)."""
        return np.array(self._snaps[0].hier.dl0)

    @property
    def domain(self):
        """Physical size of the domain, float64 (dim,)."""
        return np.array(self._snaps[0].hier.domain)

    @property
    def periodic(self):
        return np.array(self._snaps[0].hier.periodic)

    @property
    def interp_order(self):
        return None if self.meta is None else self.meta.interp_order

    def _label(self):
        return f"{type(self).__name__} {self.quantity or ''}".strip()

    # ------------------------------------------------------------------ components
    def _resolve_component(self, key):
        """Component name for ``key`` ("x", "Bx", "Pxx", "protons_Pxx"...), or None."""
        if not isinstance(key, str):
            return None
        if key in self.names:
            return key
        for name in sorted(self.names, key=len, reverse=True):
            if len(key) > len(name) and key.endswith(name):
                prefix = key[: -len(name)]
                if _is_text_prefix(prefix):
                    return name
        return None

    def _component(self, name):
        if self.names == SCALAR_NAMES:
            return self
        snaps = [
            Snapshot(s.hier, tuple(component_level(ld, name) for ld in s.levels))
            for s in self._snaps
        ]
        quantity = f"{self.quantity}_{name}" if self.quantity else name
        field = make_field(snaps, meta=self.meta, quantity=quantity, stencil=self.stencil,
                           sim=self.sim)
        field._time_dtype = self._time_dtype
        field._component_alias = name  # master's key for B.x is "x" (compat.py)
        if self._origin is not None and self._origin["component"] is None:
            field._origin = dict(self._origin, component=name)
        return field

    def __getitem__(self, key):
        """Component as a ScalarField sharing buffers (no copy, no read).

        Accepts the component name ("x", "xx") or that name after a text prefix ("Bx",
        "Pxx", "protons_Pxx"). KeyError listing the names otherwise (``B.x["x"]`` raises,
        as on master; master's per-patch keys of ``B.x`` are answered by the views).
        """
        name = self._resolve_component(key)
        if name is None:
            raise KeyError(
                f"{self._label()} has no component {key!r}; components: {self.names} "
                f"(a text prefix is accepted, e.g. 'B{self.names[0]}')"
            )
        return self._component(name)

    def __getattr__(self, attr):
        # called only when normal lookup fails: B.x, P.xx
        if attr.startswith("_") or "_snaps" not in self.__dict__:
            raise AttributeError(attr)
        if attr in self.names and self.names != SCALAR_NAMES:
            return self._component(attr)
        raise AttributeError(
            f"{type(self).__name__!r} object has no attribute {attr!r}; "
            f"components: {self.names}"
        )

    # ------------------------------------------------------------------ memory
    def nbytes(self):
        """Logical size in bytes of all times, levels and components, loaded or not."""
        return int(sum(ld.nbytes() for s in self._snaps for ld in s.levels))

    @property
    def loaded_nbytes(self):
        """Bytes of the patches currently in memory."""
        return int(sum(ld.loaded_nbytes() for s in self._snaps for ld in s.levels))

    @property
    def loaded_fraction(self):
        """``loaded_nbytes / nbytes()`` in [0, 1]."""
        total = self.nbytes()
        return self.loaded_nbytes / total if total else 1.0

    def _slots_for(self, ld, box):
        return None if box is None else ld.geom.patches_in(box)

    def load(self, box=None, force=False):
        """Read now what is not in memory yet; returns ``self``.

        box   : None (everything) or ``dict(x=(x0, x1), ...)`` / PhysBox: only the patches
                whose cells meet the box (PhysBox rule), at every level and time
        force : go ahead when the memory check says it does not fit
        One memory check for the whole request first (PhareseeMemoryError), then one
        ``ensure`` per level and time.
        """
        pbox = PhysBox.parse(box, self.ndim)
        todo = [(ld, self._slots_for(ld, pbox)) for s in self._snaps for ld in s.levels]
        needed = sum(_missing_bytes(ld, slots) for ld, slots in todo)
        if needed:
            plan.default_budget().check(needed, f"loading {self._label()}", force=force)
        for ld, slots in todo:
            if slots is None or len(slots):
                ld.ensure(slots, force=True)
        return self

    def unload(self):
        """Forget the data read from the file (next use reads again). In-memory data stays."""
        for s in self._snaps:
            for ld in s.levels:
                ld.unload()

    def copy(self):
        """Writeable in-memory copy (reads everything first, then copies every buffer)."""
        self.load()
        snaps = []
        for s in self._snaps:
            levels = tuple(
                LevelData.from_arrays(ld.geom, ld.time,
                                      {n: np.array(ld.buffer(n), copy=True) for n in ld.names})
                for ld in s.levels
            )
            snaps.append(Snapshot(s.hier, levels))
        field = self._like(snaps)
        field._component_alias = self._component_alias
        return field

    def __copy__(self):
        return self._sharing(self._snaps)

    def __deepcopy__(self, memo):
        return self.copy()

    def __reduce__(self):
        if self._origin is not None:
            o = self._origin
            return (_rebuild_from_spec, (o["spec"], o["kwargs"], o["component"]))
        states = [(s.hier.extent, tuple(s.levels)) for s in self._snaps]
        attrs = dict(meta=self.meta, quantity=self.quantity, stencil=self.stencil,
                     sim=self.sim)
        return (_rebuild_from_levels, (type(self).__name__, states, attrs))

    # ------------------------------------------------------------------ repr
    def __repr__(self):
        lines = [f"{type(self).__name__} {self.quantity or '(computed)'}  "
                 f"components {self.names}  {self.dtype}"]
        t = self._times
        if len(t) == 1:
            lines.append(f"  time {t[0]:g}")
        else:
            shown = ", ".join(f"{v:g}" for v in t[:6]) + (", ..." if len(t) > 6 else "")
            lines.append(f"  {len(t)} times [{shown}] (levels below at t={t[0]:g})")
        for ld in self._snaps[0].levels:
            g = ld.geom
            lines.append(f"  level {g.level}: {_fmt_int(g.npatch)} patches / "
                         f"{_fmt_int(g.nnodes)} nodes")
        ext = "whole domain" if self.extent is None else repr(self.extent)
        lines.append(f"  extent {ext}")
        lines.append(f"  domain {_fmt_vec(self.domain)}, dl {_fmt_vec(self.dl)}, periodic "
                     f"{tuple(bool(p) for p in self.periodic)}")
        total = self.nbytes()
        lines.append(f"  loaded {100.0 * self.loaded_fraction:.1f} % "
                     f"({_fmt_bytes(self.loaded_nbytes)} of {_fmt_bytes(total)})")
        file_backed = any(ld.backend is not None for s in self._snaps for ld in s.levels)
        if self.meta is not None:
            src = f"{self.meta.path} [{self.meta.format}]"
            lines.append(f"  file {src}" if file_backed else f"  in memory, from {src}")
            if self.meta.writer_signature:
                lines.append(f"  writer {self.meta.writer_signature}")
        else:
            lines.append("  in memory")
        if self.stencil is not None:
            lines.append(f"  stencil {self.stencil}")
        return "\n".join(lines)

    # ------------------------------------------------------------------ delegations
    def finest(self, time=None, **kw):
        """``views.finest.finest`` of the snapshot at ``time`` -> FinestGrid.

        Keywords: box, max_level, fill ("cubic" | "nearest"), dtype, force, names. Reads only
        the patches under the box. Needs one time (``time=`` or a single-time field).
        """
        from .views.finest import finest

        return finest(self.snapshot(time), **kw)

    def plot(self, **kwargs):
        """Picture of the field (``views.plot.plot``); returns ``(fig, ax)``.

        Same keywords as the legacy ``plot`` (qty, time, levels, ax, cmap, vmin, vmax,
        plot_patches, xlim, ylim, filename, ...) plus robust, symmetric, norm, max_level,
        fill. 2D: one image of ``finest()`` over ``xlim``/``ylim``; reads only the patches
        under the view. Colour limits default to the min/max of what is shown.
        """
        from .views.plot import plot

        return plot(self, **kwargs)

    def __call__(self, x=None, y=None, z=None, points=None, time=None):
        """Values at points or along a line (``views.probe.probe``), physical coordinates.

        ``B.z(y=800.)`` in 2D -> ``(x_coords, values)``; ``f(x=xs, y=ys)`` / ``f(points=P)``
        -> values (dict per component for vectors). Reads only the patches holding the points.
        """
        from .views.probe import probe

        return probe(self.snapshot(time), x=x, y=y, z=z, points=points)

    def _reduce(self, what, time, per_level):
        from .ops import reduce

        snap = self.snapshot(time)
        if per_level:
            return reduce.per_level(snap, what)
        return {"mean": reduce.mean, "min": reduce.minimum, "max": reduce.maximum}[what](snap)

    def mean(self, time=None, per_level=False):
        """Area/volume-weighted mean over the extent with the best data (float64; per
        component array for vectors). ``per_level=True``: master's per-level means (dict)."""
        return self._reduce("mean", time, per_level)

    def integral(self, time=None):
        """Integral over the extent with the best data (coarse under fine excluded), float64."""
        from .ops import reduce

        return reduce.integral(self.snapshot(time))

    def min(self, time=None, per_level=False):
        """Smallest value over the best data inside the extent (``per_level``: dict per level)."""
        return self._reduce("min", time, per_level)

    def max(self, time=None, per_level=False):
        """Largest value over the best data inside the extent (``per_level``: dict per level)."""
        return self._reduce("max", time, per_level)

    def check_shared_nodes(self, time=None):
        """Largest |difference| between copies of nodes shared by same-level patches.

        A diagnostic (reads every patch). Single-rank PHARE output has bit-identical shared
        nodes (0 everywhere); multi-rank output can differ at the 1e-4 relative level
        between copies held by patches of different MPI ranks (probes and ``finest()`` then
        show one copy). Periodic images count as shared. Returns
        ``{level: float}`` (max over components; 0.0 when a level has no shared node).
        """
        # DESIGN-NOTE: one number per level (max over components), not per component.
        snap = self.snapshot(time)
        return {ld.level: _shared_node_spread(ld) for ld in snap.levels}


def _is_text_prefix(prefix):
    """True when ``prefix`` is non-numeric text ("B", "P", "protons_P")."""
    if not prefix or not any(c.isalpha() for c in prefix):
        return False
    try:
        float(prefix)  # "inf", "nan" are numbers for float()
    except ValueError:
        return True
    return False


def _missing_bytes(ld, slots):
    """Bytes ``ld.ensure(slots)`` would read (0 for in-memory or loaded patches). No I/O."""
    return int(ld.read_nbytes(slots))


def _shared_node_spread(ld):
    """Max |difference| between copies of shared nodes of one level (all components)."""
    g = ld.geom
    if g.npatch == 0:
        return 0.0
    nlo, nhi = g.lower, g.upper + 1  # closed node boxes
    shifts = periodic_shifts(g.domain_cells, g.periodic)
    pairs = BinIndex(nlo, nhi).pairs(nlo, nhi, shifts)
    keep = (pairs.a < pairs.b) | ((pairs.a == pairs.b) & (pairs.shift != 0))
    a, b = pairs.a[keep], pairs.b[keep]
    if len(a) == 0:
        return 0.0
    s = shifts[pairs.shift[keep]]
    lo, hi = pairs.lo[keep], pairs.hi[keep]  # in b's coordinates
    size = hi - lo + 1
    ia = box_flat_indices(g.offset[a], g.shape[a], lo - s - nlo[a], size)
    ib = box_flat_indices(g.offset[b], g.shape[b], lo - nlo[b], size)
    worst = 0.0
    for f in ld.buffers().values():
        if len(ia):
            worst = max(worst, float(np.max(np.abs(f[ia].astype(np.float64)
                                                   - f[ib].astype(np.float64)))))
    return worst


def _snapshot_from_backend(backend, meta, t, box, cache, names, dtype, derived,
                           selection=None):
    """Snapshot of ``backend`` at file time ``t`` (no field bytes read).

    Patches: those of ``box`` (PhysBox rule), or of ``selection`` (MasterSelection) when given.
    """
    dl0 = np.asarray(meta.dl0, dtype=np.float64)
    cells0 = np.asarray(meta.domain_cells0, dtype=np.int64)
    geoms, lds = [], []
    for li in sorted(backend.index(t), key=lambda li: li.level):
        if li.level != len(geoms):  # a missing level number ends the hierarchy
            break
        r = 2**li.level
        dl, cells = dl0 / r, cells0 * r
        geom = cache.get(li.level, li.lower, li.upper, dl, cells, meta.periodic)
        file_ids = None
        if selection is not None:
            file_ids = np.flatnonzero(selection.selects(li.lower, li.upper, dl))
            if len(file_ids) == 0:
                break
            if len(file_ids) == geom.npatch:
                file_ids = None
            else:
                geom = cache.get(li.level, li.lower[file_ids], li.upper[file_ids], dl, cells,
                                 meta.periodic)
        elif box is not None:
            slots = geom.patches_in(box)
            if len(slots) == 0:
                break
            if len(slots) < geom.npatch:
                file_ids = np.sort(geom.file_index[slots])
                geom = cache.get(li.level, li.lower[file_ids], li.upper[file_ids], dl, cells,
                                 meta.periodic)
        geoms.append(geom)
        lds.append(LevelData.from_backend(backend, geom, t, names=names, dtype=dtype,
                                          file_ids=file_ids, derived=derived))
    if not geoms:
        raise ValueError(f"{meta.path}: no level 0 at time {t}")
    return Snapshot(HierGeometry(geoms, box), tuple(lds))


# --------------------------------------------------------------------------- public classes


class ScalarField(AMRField):
    """A one-component field (``names == ("value",)``): density, a component of B, ..."""

    _expected = (lambda names: len(names) == 1, "exactly one component")


class VectorField(AMRField):
    """A three-component field ``("x", "y", "z")``: B, E, bulk velocity, fluxes, ...

    Components: ``B.x``, ``B["x"]``, ``B["Bx"]`` (ScalarField sharing buffers).
    """

    _expected = (lambda names: tuple(names) == VECTOR_NAMES, "the components ('x', 'y', 'z')")


class TensorField(AMRField):
    """A 6-component (symmetric) or 9-component tensor field: ``P["xx"]``, ``P["Pxy"]``..."""

    _expected = (lambda names: len(names) in (6, 9), "6 (symmetric) or 9 components")


def make_field(snaps, **kw):
    """Field of the class matching the component names: 1 -> ScalarField,
    ("x", "y", "z") -> VectorField, 6 or 9 -> TensorField, else AMRField.

    Keywords: meta, quantity, stencil, sim (see AMRField).
    """
    snaps = tuple(snaps)
    if not snaps:
        raise ValueError("make_field needs at least one Snapshot")
    names = tuple(snaps[0].names)
    for cls in (ScalarField, VectorField, TensorField):
        if cls._expected[0](names):
            return cls(snaps, **kw)
    return AMRField(snaps, **kw)


__all__ = ["AMRField", "MasterSelection", "ScalarField", "TensorField", "VectorField",
           "make_field"]
