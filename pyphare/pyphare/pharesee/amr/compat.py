"""Master's per-patch hierarchy API on top of the flat per-level fields (DESIGN 8.1, 8.2).

Scripts written for master's ``PatchHierarchy`` walk levels and patches::

    for ilvl, lvl in B.levels().items():
        for patch in lvl.patches:
            pd = patch.patch_datas["x"]
            pd.dataset[ng:-ng, ng:-ng], pd.x, patch.box, patch.layout ...

``MasterAPI`` (a mixin of ``field.AMRField``) answers that API with lazy views:

- ``levels()`` builds one small ``LevelView`` per level; no patch object exists yet.
- ``LevelView.patches`` is a lazy sequence in FILE order (vtkhdf: ``AMRBox`` row order;
  phareh5: groups ``p<rank>#<localid>`` sorted by ``(rank, localid)``): ``patches[i]`` is
  the same box in both formats. A ``PatchView`` is built when it is indexed.
- ``patch.patch_datas[name]`` is a ``FieldDataView`` (a master ``FieldData``) whose
  ``dataset`` is built on first access: a READ-ONLY COPY shaped like master's
  ``all_primal=True`` output, the patch's primal nodes ``lower .. upper + 1`` surrounded by
  a NaN ring of ``ghosts_nbr`` nodes on each side of every axis, indexed ``[x, y, z]``.
  The first ``dataset`` touched on a level reads the WHOLE level in one backend call when
  it fits the memory budget (below the warning fraction of ``plan.MemoryBudget``), so that
  master's loop over every patch costs one read per level instead of one file open per
  patch; when the level does not fit, only that patch is read.

Ghost width (``ghosts_nbr``): phareh5 files = the dataset's ``ghosts`` attribute, as master
(fromh5.py:117-120); otherwise (vtkhdf, in-memory results) ``[2, 4, 4][interp_order - 1]``
for hybrid fields (gridlayout.py:225) and 2 when the field has no ``interp_order`` (MHD).

Documented differences with master (the differential table of DESIGN 11.4):

- ``time_hier`` keys are floats, not ``"%.10f"`` strings; master's formatted string is
  still accepted as a lookup key everywhere a time is taken.
- a multi-time field never picks "the first time" silently: ``levels()``, ``level()``,
  ``levelNbr()``... without ``time=`` raise ``ValueError`` listing the times.
- ghost nodes are NaN for every quantity (master had real ghost values for moments);
  ``FieldDataView.compare`` (so ``hierarchy_compare``) treats NaN == NaN at the same place.
- ``patch.id`` on vtkhdf files is ``"vtk#<file index>"`` (no MPI rank in the file: master's
  ``int(patch.id.split("#")[0][1:])`` raises instead of answering rank 0); in-memory results
  use ``"mem#<index>"``.
- a component field (``B.x``) answers master's keys: ``patch_datas`` and ``quantities()``
  list ``"x"`` (master's sub-hierarchy key) and ``"value"`` is accepted too.
- ``dataset`` keeps the field's dtype (master returned float64 for B and E).
- ``data_files`` and ``selection_box`` are empty dicts.

Views are cached weakly per (field, time): they die with the field unless the caller
keeps them. ``pyphare.pharesee.hierarchy`` and ``pyphare.core`` are imported inside
functions only (master's hierarchy package imports this package, not the reverse).
"""

from __future__ import annotations

import functools
import weakref
from collections.abc import Mapping, Sequence

import numpy as np

# hybrid field ghost width per interpolation order, as master's GridLayout default
# (pyphare/core/gridlayout.py:225) and the C++ layout
_HYBRID_FIELD_GHOSTS = (2, 4, 4)
# width used when the field has no interpolation order (MHD vtkhdf files)
_NO_INTERP_GHOSTS = 2

# field -> {time: (snapshot levels tuple, {ilvl: LevelView})}; weak on the field so that
# views die with it
_VIEW_CACHE = weakref.WeakKeyDictionary()


def default_ghost_width(interp_order):
    """Ghost width master shows for a field without a per-dataset ``ghosts`` attribute.

    interp_order : 1, 2, 3 (hybrid) or None (MHD). Returns ``[2, 4, 4][interp_order - 1]``,
    or 2 for None. Raises ValueError for another interpolation order.
    """
    if interp_order is None or int(interp_order) == 0:
        return _NO_INTERP_GHOSTS
    if int(interp_order) not in (1, 2, 3):
        raise ValueError(f"interp_order={interp_order}: expected 1, 2, 3 or None")
    return _HYBRID_FIELD_GHOSTS[int(interp_order) - 1]


class _FieldInfo:
    """What the views need from the field, copied so that views never hold the field."""

    __slots__ = ("alias", "interp_order", "keys", "meta", "names", "ndim")

    def __init__(self, field):
        self.names = tuple(field.names)
        self.meta = field.meta
        self.interp_order = field.interp_order
        self.ndim = int(field.ndim)
        # a component field (B.x) is shown with its parent's component name, as master's
        # sub-hierarchies (hierarchy.py:72-86): keys ("x",), storage name "value"
        self.alias = _component_alias(field)
        self.keys = (self.alias,) if self.alias is not None else self.names

    def storage_name(self, key):
        """Component name in the store for a ``patch_datas`` key, or None if unknown."""
        if key in self.names:
            return key
        if self.alias is not None and key == self.alias:
            return self.names[0]
        return None

    def label(self, name):
        """The key master shows for the stored component ``name``."""
        return self.alias if self.alias is not None else name


def _component_alias(field):
    """Parent component name of a component field (``"x"`` for ``B.x``), else None."""
    alias = getattr(field, "_component_alias", None)
    if alias is None or tuple(field.names) != ("value",):
        return None
    return alias


_PICKLE_VIEW_MSG = (
    "{what} cannot be pickled: it is a lazy view of a field. Pickle the field itself "
    "(pickle.dumps(field), small for a file-backed field), or copy.deepcopy(view) for a "
    "master object holding the data"
)


class _LevelContext:
    """Shared state of the views of one level at one time (one LevelData).

    Tables are built on first use and are vectorised over the level's patches; the
    patch views themselves are kept in a WeakValueDictionary (file-order index -> view).
    """

    def __init__(self, info, ld):
        self.info = info
        self.ld = ld
        self.geom = ld.geom
        self.level = int(ld.geom.level)
        self.time = float(ld.time)
        self._boxes = None
        self._payload = None
        self._payload_done = False
        self._ghosts = {}
        self._views = weakref.WeakValueDictionary()

    @property
    def npatch(self):
        return int(self.geom.npatch)

    def boxes(self):
        """(lower, upper) int64 ``(npatch, dim)`` in FILE order."""
        if self._boxes is None:
            self._boxes = self.geom.file_order_boxes()
        return self._boxes

    def backend_ids(self):
        """Backend file-order patch id of each file-order patch of this (sub-)level."""
        if self.ld.file_ids is not None:  # a sub-level: file_ids[i] is the backend's id
            return np.asarray(self.ld.file_ids, dtype=np.int64)
        return np.arange(self.npatch, dtype=np.int64)

    def payload(self):
        """The backend's index payload for this level (H5LevelTable on phareh5), or None.

        Reads the backend's in-memory index (already built when the field was opened).
        """
        if not self._payload_done:
            backend = self.ld.backend
            if backend is not None:
                for li in backend.index(self.time):
                    if li.level == self.level:
                        self._payload = li.payload
                        break
            self._payload_done = True
        return self._payload

    def ghosts(self, name):
        """Ghost width of component ``name`` for every file-order patch, int64 ``(npatch,)``."""
        got = self._ghosts.get(name)
        if got is not None:
            return got
        payload = self.payload()
        table = getattr(payload, "ghosts", None)
        meta = self.info.meta
        if table is not None and meta is not None:
            names = tuple(meta.names)
            # DESIGN-NOTE: a component that is not one of the file's (a backend-derived
            # quantity) takes the ghost width of the file's first dataset.
            c = names.index(name) if name in names else 0
            got = np.asarray(table)[c][self.backend_ids()].astype(np.int64)
        else:
            got = np.full(self.npatch, default_ghost_width(self.info.interp_order), np.int64)
        self._ghosts[name] = got
        return got

    def patch(self, i):
        """PatchView of file-order patch ``i`` (negative indexes allowed)."""
        n = self.npatch
        if not isinstance(i, (int, np.integer)):
            raise TypeError(f"patch index must be an int, got {type(i).__name__}")
        i = int(i)
        if i < 0:
            i += n
        if not 0 <= i < n:
            raise IndexError(f"patch index {i} out of range: level {self.level} has {n} patches")
        view = self._views.get(i)
        if view is None:
            view = PatchView(self, i)
            self._views[i] = view
        return view


class _PatchSequence(Sequence):
    """Lazy, read-only sequence of the PatchViews of a level, in FILE order.

    Supports ``len``, ``p[i]`` (negative too), ``p[a:b]`` (a list), iteration; a view is
    built only when its index is accessed.
    """

    __slots__ = ("_ctx",)

    def __init__(self, ctx):
        self._ctx = ctx

    def __len__(self):
        return self._ctx.npatch

    def __getitem__(self, idx):
        if isinstance(idx, slice):
            return [self._ctx.patch(i) for i in range(*idx.indices(self._ctx.npatch))]
        return self._ctx.patch(idx)

    def __iter__(self):
        for i in range(self._ctx.npatch):
            yield self._ctx.patch(i)

    def __repr__(self):
        return f"<{self._ctx.npatch} patches of level {self._ctx.level} (lazy)>"


class LevelView:
    """Master's ``PatchLevel`` for one level of a field at one time (lazy).

    level_number : int
    patches      : lazy sequence of PatchView in FILE order (no patch object until accessed)
    cell_width   : (ndim,) float64 mesh size of this level
    Reads nothing from disk.
    """

    __slots__ = ("__weakref__", "_ctx", "level_number", "patches")

    def __init__(self, ctx):
        self._ctx = ctx
        self.level_number = ctx.level
        self.patches = _PatchSequence(ctx)

    def __iter__(self):
        return iter(self.patches)

    def __getitem__(self, idx):
        if isinstance(idx, (int, np.integer)):
            return self.patches[idx]
        raise IndexError(f"PatchLevel::__getitem__ unhandled input type: {type(idx)}")

    @property
    def cell_width(self):
        """Mesh size of this level, float64 ``(ndim,)``."""
        return np.array(self._ctx.geom.dl, dtype=np.float64)

    def level_range(self):
        """``(xmin, xmax)`` of the node coordinates along x, ghost nodes included (master)."""
        ctx = self._ctx
        lower, upper = ctx.boxes()
        g = ctx.ghosts(ctx.info.names[0])
        dx = float(ctx.geom.dl[0])
        return float(((lower[:, 0] - g) * dx).min()), float(((upper[:, 0] + 1 + g) * dx).max())

    def __deepcopy__(self, memo):
        """Master ``PatchLevel`` of master ``Patch`` copies (reads the whole level)."""
        from pyphare.pharesee.hierarchy.patchlevel import PatchLevel

        return PatchLevel(self.level_number, [p.copy() for p in self.patches])

    def __copy__(self):
        return LevelView(self._ctx)

    def __reduce__(self):
        raise TypeError(_PICKLE_VIEW_MSG.format(what="a LevelView (field.level(i))"))

    def __repr__(self):
        return f"LevelView(level={self.level_number}, {self._ctx.npatch} patches)"


class PatchView:
    """Master's ``Patch`` for one patch, built when the patch is indexed; reads nothing.

    box         : pyphare Box, inclusive cell indexes of this level
    origin      : float64 ``(ndim,)`` physical position of node ``lower``
    dl          : float64 ``(ndim,)`` mesh size of the level
    id          : phareh5: the HDF5 group name ``p<rank>#<localid>``; other formats:
                  ``"p0#<i>"`` with ``i`` the file-order index
    attrs       : dict with ``lower``, ``upper``, ``nbrCells``, ``origin`` (+ ``mpi_rank`` on
                  phareh5, parsed from the group name)
    layout      : GridLayout of the first component (master's ``patch.layout``)
    patch_datas : dict-like ``{name: FieldDataView}``, values built on access
    """

    __slots__ = ("__weakref__", "_box", "_ctx", "_i", "_pdatas")

    created = 0  # instances built since import (tests check laziness with it)

    def __init__(self, ctx, i):
        self._ctx = ctx
        self._i = int(i)
        self._box = None
        self._pdatas = None
        PatchView.created += 1

    @property
    def box(self):
        if self._box is None:
            from pyphare.core.box import Box

            lower, upper = self._ctx.boxes()
            self._box = Box(lower[self._i], upper[self._i])
        return self._box

    @property
    def origin(self):
        lower, _ = self._ctx.boxes()
        return lower[self._i] * np.asarray(self._ctx.geom.dl, dtype=np.float64)

    @property
    def dl(self):
        return np.array(self._ctx.geom.dl, dtype=np.float64)

    @property
    def id(self):
        payload = self._ctx.payload()
        if hasattr(payload, "patch_name"):
            return payload.patch_name(int(self._ctx.backend_ids()[self._i]))
        # DESIGN-NOTE: vtkhdf has no rank nor patch name. The id keeps master's "<x>#<i>"
        # shape but must not parse as a rank (master's int(id.split("#")[0][1:]) would count
        # every patch on rank 0): "vtk#<file index>" for a vtkhdf file, "mem#<index>" for
        # in-memory results.
        ld = self._ctx.ld
        if ld.backend is not None and getattr(ld.backend.meta(), "format", None) == "vtkhdf":
            return f"vtk#{int(self._ctx.backend_ids()[self._i])}"
        return f"mem#{self._i}"

    @property
    def attrs(self):
        lower, upper = self._ctx.boxes()
        lo, up = lower[self._i].copy(), upper[self._i].copy()
        attrs = {
            "lower": lo,
            "upper": up,
            "nbrCells": up - lo + 1,
            "origin": self.origin,
        }
        payload = self._ctx.payload()
        rank = getattr(payload, "rank", None)
        if rank is not None:
            attrs["mpi_rank"] = np.uint64(rank[int(self._ctx.backend_ids()[self._i])])
        return attrs

    @property
    def patch_datas(self):
        if self._pdatas is None:
            self._pdatas = _PatchDatas(self._ctx, self._i)
        return self._pdatas

    @property
    def layout(self):
        return self.patch_datas[self._ctx.info.names[0]].layout

    def __getitem__(self, key):
        if key in self.patch_datas:
            return self.patch_datas[key]
        raise KeyError(f"No patchdata for key: {key} in {list(self.patch_datas)}")

    def __call__(self, qty, **kwargs):
        """Master's 2D cut ``patch(qty, x=...)`` / ``patch(qty, y=...)`` (patch.py:51-68)."""
        from pyphare.pharesee.hierarchy.patch import Patch

        return Patch.__call__(self, qty, **kwargs)

    def copy(self):
        """Master ``Patch`` holding writeable copies of every component (reads the patch)."""
        import copy

        from pyphare.pharesee.hierarchy.patch import Patch

        pdatas = {k: copy.deepcopy(v) for k, v in self.patch_datas.items()}
        return Patch(pdatas, self.id, box=self.box, attrs=self.attrs)

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memo):
        # master code deep-copies patches to mutate them (geometry.shift_patch): give it a
        # master Patch; the view itself stays read-only
        return self.copy()

    def __reduce__(self):
        raise TypeError(_PICKLE_VIEW_MSG.format(what="a PatchView (level.patches[i])"))

    def __str__(self):
        return f"Patch: box( {self.box}), id({self.id})"

    def __repr__(self):
        return self.__str__()


class _PatchDatas(Mapping):
    """``patch.patch_datas``: component name -> FieldDataView, each built on first access."""

    __slots__ = ("_built", "_ctx", "_i")

    def __init__(self, ctx, i):
        self._ctx = ctx
        self._i = i
        self._built = {}

    def __getitem__(self, key):
        info = self._ctx.info
        name = info.storage_name(key) if isinstance(key, str) else None
        if name is None:
            raise KeyError(key)
        pd = self._built.get(name)
        if pd is None:
            pd = _make_field_data_view(self._ctx, self._i, name)
            self._built[name] = pd
        return pd

    def __contains__(self, key):
        return isinstance(key, str) and self._ctx.info.storage_name(key) is not None

    def __iter__(self):
        return iter(self._ctx.info.keys)

    def __len__(self):
        return len(self._ctx.info.keys)

    def __repr__(self):
        return f"<patch_datas {list(self._ctx.info.keys)} (lazy)>"


@functools.cache
def _field_data_view_class():
    """FieldDataView, created on first use (it subclasses master's FieldData, imported late)."""
    from pyphare.pharesee.hierarchy.patchdata import FieldData

    class FieldDataView(FieldData):
        """Master ``FieldData`` of one component of one patch, data read on demand.

        Same attributes as master's ``all_primal=True`` FieldData: ``layout`` (a real
        GridLayout), ``box``, ``origin``, ``dl``, ``ndim``, ``field_name``/``name``,
        ``centerings`` (all "primal"), ``ghosts_nbr``, ``ghost_box``, ``size``,
        ``x/y/z`` (node coordinates INCLUDING the ghost ring); ``select``,
        ``__getitem__(box)``, ``meshgrid``, ``compare`` are master's.
        ``dataset``: read-only ndarray copy, shape ``size``, indexed ``[x, y, z]``,
        primal nodes inside, NaN in the ghost ring; built on first access (reads this
        patch only) and kept. Assigning ``dataset`` raises AttributeError, writing into
        it raises ValueError; ``copy.deepcopy`` gives a writeable master FieldData.
        """

        created = 0

        def __init__(self, ctx, i, name):
            from pyphare.core import box as boxm
            from pyphare.core.gridlayout import GridLayout

            lower, upper = ctx.boxes()
            ndim = ctx.info.ndim
            dl = np.array(ctx.geom.dl, dtype=np.float64)
            ng = int(ctx.ghosts(name)[i])
            ghosts = np.full(ndim, ng, dtype=int)
            layout = GridLayout(
                boxm.Box(lower[i], upper[i]), lower[i] * dl, dl,
                0 if ctx.info.interp_order is None else int(ctx.info.interp_order),
                [ng] * ndim,
            )
            # the attributes of PatchData.__init__ and FieldData.__init__
            self.quantity = "field"
            self.box = layout.box
            self.origin = layout.origin
            self.layout = layout
            self._x = None
            self._y = None
            self._z = None
            self.dl = np.asarray(layout.dl)
            label = ctx.info.label(name)  # "x" for B.x, as master's sub-hierarchy
            self.field_name = label
            self.name = label
            self.ndim = ndim
            self.centerings = ["primal"] * ndim
            self.ghosts_nbr = ghosts
            self.ghost_box = boxm.grow(layout.box, ghosts)
            self.size = np.asarray(self.ghost_box.shape) + 1
            self.offset = np.zeros(ndim)
            self._src = (ctx.ld, int(ctx.geom.store_index[i]), name)
            self._data = None
            type(self).created += 1

        @property
        def dataset(self):
            if self._data is None:
                ld, slot, name = self._src
                _load_level_if_it_fits(ld, slot)
                self._data = _ringed_copy(ld.view(slot, name), int(self.ghosts_nbr[0]))
            return self._data

        @dataset.setter
        def dataset(self, value):
            raise AttributeError(
                "the dataset of a file-backed patch view is read-only; use "
                "copy.deepcopy(patch_data) or field.copy() for writeable data"
            )

        def __deepcopy__(self, memo):
            import copy

            return FieldData(
                copy.deepcopy(self.layout, memo), self.field_name,
                np.array(self.dataset), centering=list(self.centerings),
                ghosts_nbr=np.array(self.ghosts_nbr),
            )

        def __reduce__(self):
            return self.__deepcopy__({}).__reduce__()

        def compare(self, that, atol=1e-16):
            """Master's ``FieldData.compare`` with NaN == NaN at the same place.

            True when the names match and the datasets are equal within ``atol`` (master's
            float32/float64 tolerance rule, ``phare_utilities.any_fp_tol``), a NaN in one
            matching a NaN at the same index in the other: the NaN ghost ring of two views
            of the same data compares equal (``hierarchy_compare(Ni, Ni_again)`` is True).
            Reads both datasets.
            """
            from pyphare.core import phare_utilities as phut

            a, b = self.dataset[:], that.dataset[:]
            if self.field_name != that.field_name or np.shape(a) != np.shape(b):
                return False
            return bool(np.allclose(a, b, equal_nan=True,
                                    **phut.any_fp_tol(a, b, atol, 0, None)))

        def __eq__(self, that):
            return self.compare(that)

        def __ne__(self, that):
            return not self.compare(that)

        __hash__ = None  # as master's FieldData (it defines __eq__)

    return FieldDataView


def _load_level_if_it_fits(ld, slot):
    """Read the whole level of ``ld`` at once when patch ``slot`` is missing and it fits.

    Master's per-patch loops touch every patch of a level in turn; reading them one by one
    costs one file open and one small read each (~1.2 ms per patch on phareh5). When the
    missing bytes of the whole level are below the warning fraction of the memory budget
    (``plan.WARN_FRACTION`` of what is available), the level is read in one backend call;
    otherwise nothing is done here and the caller reads the single patch. Never raises for
    memory: the fallback is the single-patch read.
    """
    from . import plan

    if ld.backend is None or ld.loaded[slot]:
        return
    needed = plan.ReadPlan().add_read(ld).total
    budget = ld.budget if ld.budget is not None else plan.default_budget()
    if needed < plan.WARN_FRACTION * budget.available():
        ld.ensure()


def _make_field_data_view(ctx, i, name):
    return _field_data_view_class()(ctx, i, name)


def field_data_view_class():
    """The FieldDataView class (a subclass of master's FieldData)."""
    return _field_data_view_class()


def _ringed_copy(nodes, ng):
    """Read-only copy of ``nodes`` (patch primal nodes, ``[x, y, z]``) inside a NaN ring.

    The result has ``ng`` extra nodes on both sides of every axis, filled with NaN. Float
    dtypes are kept; other dtypes become float64 (NaN needs a float).
    """
    dtype = nodes.dtype if nodes.dtype.kind in "fc" else np.dtype(np.float64)
    out = np.full(tuple(n + 2 * ng for n in nodes.shape), np.nan, dtype=dtype)
    out[tuple(slice(ng, ng + n) for n in nodes.shape)] = nodes
    out.flags.writeable = False
    return out


class _TimeHier(Mapping):
    """``hier.time_hier``: float time -> ``{ilvl: LevelView}``, built per time on lookup.

    Keys are floats (master: ``"%.10f"`` strings); a string key such as
    ``format_timestamp(t)`` is accepted for lookups and ``in``.
    """

    __slots__ = ("_field",)

    def __init__(self, field):
        self._field = field

    def _times(self):
        return [float(t) for t in np.asarray(self._field.times()).ravel()]

    def __iter__(self):
        return iter(self._times())

    def __len__(self):
        return len(self._times())

    def __getitem__(self, key):
        try:
            return self._field.levels(_as_time(key))
        except ValueError as e:  # an unparsable key is a missing key for a mapping
            raise KeyError(key) from e

    def __repr__(self):
        return f"<time_hier, times {self._times()} (lazy)>"


def _as_time(time):
    """float from a float or master's ``"%.10f"`` time string; None stays None."""
    if time is None:
        return None
    if isinstance(time, bytes):
        time = time.decode()
    if isinstance(time, str):
        try:
            return float(time)
        except ValueError:
            raise ValueError(f"time {time!r} is not a number") from None
    return float(time)


class MasterAPI:
    """Master's ``PatchHierarchy`` API on an AMR field, as lazy views (mixin of AMRField).

    Uses only ``self.times()``, ``self.snapshot(time)``, ``self.names``, ``self.meta``,
    ``self.interp_order``, ``self.ndim``, ``self.dl``, ``self.domain``.
    ``time`` arguments take a float (matched with the backend's tolerance) or master's
    ``"%.10f"`` string. A field holding several times raises ValueError (listing them)
    where master would have used the first time silently. Nothing here reads field data
    except ``global_min``/``global_max`` and the ``dataset`` of the patch views.
    """

    refinement_ratio = 2

    # ------------------------------------------------------------------ time helpers
    def _compat_time(self, time):
        """The float time a structural call works on (``time=None`` needs a single time)."""
        time = _as_time(time)
        if time is None:
            times = [float(t) for t in np.asarray(self.times()).ravel()]
            if len(times) != 1:
                raise ValueError(
                    f"this field holds {len(times)} times {times}: pass time=... "
                    "(master's PatchHierarchy used the first one silently)"
                )
            time = times[0]
        return time

    def _default_time(self):
        """The field's only time (master's ``_default_time``); ValueError if several."""
        return self._compat_time(None)

    def _compat_levels(self, time):
        snap = self.snapshot(self._compat_time(time))
        key = float(snap.time)
        try:
            per_field = _VIEW_CACHE.setdefault(self, {})
        except TypeError:  # field not weak-referenceable: no cache
            per_field = {}
        cached = per_field.get(key)
        if cached is not None and cached[0] is snap.levels:
            return cached[1]
        info = _FieldInfo(self)
        levels = {int(ld.geom.level): LevelView(_LevelContext(info, ld)) for ld in snap.levels}
        per_field[key] = (snap.levels, levels)
        return levels

    # ------------------------------------------------------------------ structure
    def levels(self, time=None):
        """``{ilvl: LevelView}`` at ``time`` (master: ``{ilvl: PatchLevel}``). No I/O.

        time : float or ``"%.10f"`` string; None on a single-time field.
        Raises ValueError on a multi-time field without ``time``, TimeNotFoundError for a
        time that is not in the field.
        """
        return dict(self._compat_levels(time))

    def level(self, level_number, time=None):
        """LevelView of ``level_number`` at ``time`` (KeyError if the level is absent)."""
        return self._compat_levels(time)[level_number]

    def levelNbr(self, time=None):
        """Number of levels at ``time``."""
        return len(self._compat_levels(time))

    def levelNbrs(self, time=None):
        """List of the level numbers at ``time``."""
        return list(self._compat_levels(time).keys())

    def finest_level(self, time=None):
        """Largest level number at ``time``."""
        return max(self._compat_levels(time).keys())

    @property
    def time_hier(self):
        """Mapping float time -> ``{ilvl: LevelView}`` (lazy; string keys accepted)."""
        return _TimeHier(self)

    def has_time(self, time):
        """True if ``time`` (float or ``"%.10f"`` string) is one of the field's times."""
        try:
            self.snapshot(_as_time(time))
        except (KeyError, ValueError):
            return False
        return True

    def has_file(self, filename):
        """Always False: a field keeps no open file (``data_files`` is empty)."""
        return False

    def quantities(self):
        """Component names, as master's ``patch_datas`` keys (``["x", "y", "z"]``...).

        A component field (``B.x``) answers ``["x"]`` like master's sub-hierarchy.
        """
        alias = _component_alias(self)
        return [alias] if alias is not None else list(self.names)

    def is_homogeneous(self):
        """Always True: every patch holds every component."""
        return True

    def nbrPatches(self):
        """Number of patches summed over all levels and all times (master's meaning)."""
        n = 0
        for t in np.asarray(self.times()).ravel():
            n += sum(int(ld.geom.npatch) for ld in self.snapshot(float(t)).levels)
        return n

    def update(self):
        """Does nothing (master rebuilt per-component sub-hierarchies here)."""

    def add_time(self, *args, **kwargs):
        """Not supported: AMR fields are immutable; read the times you need with Run.Get*."""
        raise NotImplementedError(
            "add_time: fields of pyphare.pharesee.amr are immutable; ask Run.GetX for the "
            "times you need (a list of times is accepted), or use Run(path, legacy=True)"
        )

    # ------------------------------------------------------------------ boxes
    @property
    def domain_box(self):
        """Level-0 cell Box ``[0, N0 - 1]`` of the whole domain."""
        from pyphare.core.box import Box

        cells = np.rint(np.asarray(self.domain, dtype=np.float64)
                        / np.asarray(self.dl, dtype=np.float64)).astype(np.int64)
        return Box(np.zeros_like(cells), cells - 1)

    def refined_domain_box(self, level_number):
        """Domain Box refined to ``level_number``."""
        from pyphare.core import box as boxm

        assert level_number >= 0
        return boxm.refine(self.domain_box, self.refinement_ratio**level_number)

    def level_domain_box(self, level_number):
        """Domain Box at ``level_number`` (level 0: ``domain_box``)."""
        if level_number == 0:
            return self.domain_box
        return self.refined_domain_box(level_number)

    @property
    def data_files(self):
        """Empty dict: the field opens its file per read and keeps no handle."""
        return {}

    @property
    def selection_box(self):
        """Empty dict (master: ``{time: Box}``); the field's extent replaces it."""
        # DESIGN-NOTE: a field read with box= is described by its extent (field.py), not by
        # master's selection_box Boxes; the empty dict is master's "no selection" value.
        return {}

    # ------------------------------------------------------------------ reductions
    def _compat_reduce(self, fn, qty, time):
        snap = self.snapshot(self._compat_time(time))
        alias = _component_alias(self)
        if alias is not None and qty == alias:  # B.x.global_max("x"), as master
            qty = snap.names[0]
        if qty not in snap.names:
            keys = list(snap.names) + ([alias] if alias is not None else [])
            raise KeyError(f"no component {qty!r}; this field has {keys}")
        values = [fn(ld.buffer(qty)) for ld in snap.levels if ld.geom.nnodes]
        return fn(np.asarray(values))

    def global_min(self, qty, time=None):
        """NaN-ignoring minimum of component ``qty`` over every node of every level.

        ``qty``: a component name; on a component field (``B.x``) both "x" and "value".

        Master's value (its NaN ghost ring is ignored there too). Reads the whole time.
        """
        return self._compat_reduce(np.nanmin, qty, time)

    def global_max(self, qty, time=None):
        """NaN-ignoring maximum of component ``qty`` over every node of every level."""
        return self._compat_reduce(np.nanmax, qty, time)

    # ------------------------------------------------------------------ master's plots
    def box_to_Rectangle(self, box):
        """matplotlib Rectangle of a cell Box (master)."""
        from matplotlib.patches import Rectangle

        return Rectangle(box.lower, *box.shape)

    def plot_patches(self, save=False):
        """Master's 1D patch-outline plot (it uses time 0.0, as master)."""
        from pyphare.pharesee.hierarchy.hierarchy import PatchHierarchy

        return PatchHierarchy.plot_patches(self, save=save)

    def plot_2d_patches(self, ilvl, collections, **kwargs):
        """Master's plot of Box collections on the level ``ilvl`` grid; returns the figure."""
        from pyphare.pharesee.hierarchy.hierarchy import PatchHierarchy

        return PatchHierarchy.plot_2d_patches(self, ilvl, collections, **kwargs)

    # DESIGN-NOTE: plot() belongs to the plot layer; master's per-patch plot1d/plot2d are
    # kept here, running master's code on the views (they read every patch shown).
    def plot1d(self, **kwargs):
        """Master's per-patch 1D plot (hierarchy.py:396-447) on the views."""
        from pyphare.pharesee.hierarchy.hierarchy import PatchHierarchy

        return PatchHierarchy.plot1d(self, **kwargs)

    def plot2d(self, **kwargs):
        """Master's per-patch 2D pcolormesh plot (hierarchy.py:460-541); ``(fig, ax)``."""
        from pyphare.pharesee.hierarchy.hierarchy import PatchHierarchy

        return PatchHierarchy.plot2d(self, **kwargs)

    def _pcolormesh_coords(self, pdat, box):
        # helper of master's plot2d, called on self
        from pyphare.pharesee.hierarchy.hierarchy import PatchHierarchy

        return PatchHierarchy._pcolormesh_coords(self, pdat, box)
