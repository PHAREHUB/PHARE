"""Elementwise arithmetic between AMR snapshots (DESIGN 2 "arithmetic", 7.1).

Every function takes and returns ``Snapshot`` objects (one quantity, one time: a
``HierGeometry`` plus one ``LevelData`` per level). Results are in memory
(``LevelData.from_arrays``), node centred, x fastest inside each patch, store order.

Two operands on the SAME layout (equal level keys) are combined with one numpy call per
level and per component on the whole flat buffers: no Python object per patch.

Two operands on DIFFERENT layouts (another time, another box) are first brought onto a
common layout, level by level:
  * a level whose geometry keys are equal is used as is;
  * otherwise the result level is made of the non-empty cell-box intersections of the
    operands' patches (file order of the result = order of the (a-patch, b-patch) pairs),
    and every operand is gathered onto it with flat index blocks;
  * a level present in only one operand (or whose intersection is empty) is dropped,
    together with every finer level.
Level 0 covers the whole domain in both operands, so it never loses anything. The extent of
the result is the intersection of the operands' extents.
Between two times: ``B1 - B0`` lives on the common refinement of the two layouts, so
``integral(B1 - B0)`` integrates over that layout's best data and can differ from
``integral(B1) - integral(B0)`` (each over its own best data) where the refined regions
moved; restrict both to the common layout first to compare like with like.
Operands on different grids (mesh size, number of cells or periodicity of level 0) raise
``ValueError``.

Component rules: equal names -> name by name; one operand with ``("value",)`` (a scalar
field) is broadcast over the other's names; anything else raises ``TypeError``. Python and
numpy scalars (and 0-d arrays) broadcast everywhere; arrays with ``ndim > 0`` raise
``TypeError``. dtypes follow numpy (``np.result_type`` / NEP 50): a float32 field combined
with a Python float stays float32, with a ``np.float64`` scalar it becomes float64.

Memory: the result size (and any file read) is checked with ``plan`` before allocating.
Gathers work in blocks of at most ``BLOCK_NODES`` nodes (index and value temporaries of
~32 MB each).
"""

from __future__ import annotations

import numbers

import numpy as np

from ..boxes import box_flat_indices, intersect_tables
from ..geometry import HierGeometry, default_cache
from ..plan import ReadPlan
from ..store import LevelData, Snapshot

SCALAR_NAMES = ("value",)
VECTOR_NAMES = ("x", "y", "z")

# DESIGN-NOTE: DESIGN 7.1 bounds gather blocks to 16 M nodes; the task also asks for per-node
# temporaries of ~64 MB at most, so blocks hold 4 M nodes (int64 indices = 32 MB).
BLOCK_NODES = 4 * 2**20


# --------------------------------------------------------------------------- operands


def _check_operand(x, what):
    """Return ``x`` if it is a Snapshot or a scalar; TypeError otherwise."""
    if isinstance(x, Snapshot):
        return x
    if isinstance(x, np.ndarray):
        if x.ndim > 0:
            raise TypeError(
                f"{what}: a numpy array of shape {x.shape} cannot be combined with an AMR "
                "field (its patches have no array shape). Use a scalar, another field, or "
                "turn the field into a regular array first with finest()."
            )
        return x
    if isinstance(x, (numbers.Number, np.generic)):
        return x
    raise TypeError(
        f"{what}: cannot combine an AMR field with {type(x).__name__}; "
        "use a number or another field of the same run"
    )


def _snapshots(operands):
    return [o for o in operands if isinstance(o, Snapshot)]


def result_names(snaps, what="operation"):
    """Component names of a result combining these snapshots (rules in the module doc)."""
    names = [tuple(s.names) for s in snaps]
    full = [n for n in names if n != SCALAR_NAMES]
    if not full:
        return SCALAR_NAMES
    if any(n != full[0] for n in full[1:]):
        raise TypeError(
            f"{what}: cannot combine fields with components {full[0]} and "
            f"{next(n for n in full if n != full[0])}; combine components one by one "
            "(e.g. B.x * E.y) or fields of the same kind"
        )
    return full[0]


def _component_of(ld, name):
    """Component of ``ld`` used for result component ``name`` (scalar fields broadcast)."""
    return "value" if ld.names == SCALAR_NAMES else name


# --------------------------------------------------------------------------- layouts


def check_same_grid(hiers):
    """Raise ValueError unless every HierGeometry has the same level-0 grid.

    Same grid = same dimension, level-0 mesh size, number of cells and periodicity.
    """
    h0 = hiers[0]
    for h in hiers[1:]:
        if (
            h.dim != h0.dim
            or not np.array_equal(h.dl0, h0.dl0)
            or not np.array_equal(h.domain_cells0, h0.domain_cells0)
            or not np.array_equal(h.periodic, h0.periodic)
        ):
            raise ValueError(
                "cannot combine fields on different grids: "
                f"dl {tuple(h0.dl0)} vs {tuple(h.dl0)}, cells {tuple(h0.domain_cells0)} vs "
                f"{tuple(h.domain_cells0)}, periodic {tuple(h0.periodic)} vs "
                f"{tuple(h.periodic)}; both fields must come from the same simulation grid"
            )


def common_extent(hiers):
    """Intersection of the extents (None = whole domain). ValueError when it is empty."""
    ext, seen = None, False
    for h in hiers:
        if h.extent is None:
            continue
        if not seen:
            ext, seen = h.extent, True
            continue
        ext = ext.intersect(h.extent)
        if ext is None:
            raise ValueError(
                "the fields' boxes do not overlap: nothing is left to combine; "
                "read both fields with overlapping box=..."
            )
    return ext


def intersect_levels(ga, gb):
    """LevelGeometry made of the non-empty cell-box intersections of two levels' patches.

    ``ga`` itself when the keys are equal; None when no pair of patches intersects. The
    result's file order is the order of the (a-patch, b-patch) pairs, a and b in file order.
    """
    if ga is gb or ga.key() == gb.key():
        return ga
    alo, ahi = ga.file_order_boxes()
    blo, bhi = gb.file_order_boxes()
    pairs = intersect_tables(alo, ahi, blo, bhi)
    if len(pairs.a) == 0:
        return None
    return default_cache().get(ga.level, pairs.lo, pairs.hi, ga.dl, ga.domain_cells,
                               ga.periodic)


def common_hier(hiers):
    """The HierGeometry on which several snapshots can be combined (rules in module doc).

    Returns the first one itself when every level is shared and the extent is unchanged.
    Raises ValueError for different grids or extents that do not overlap.
    """
    check_same_grid(hiers)
    extent = common_extent(hiers)
    h0 = hiers[0]
    levels = []
    for g0 in h0.levels:
        ilvl = g0.level
        if not all(h.has_level(ilvl) for h in hiers[1:]):
            break
        g = g0
        for h in hiers[1:]:
            g = intersect_levels(g, h.level(ilvl))
            if g is None:
                break
        if g is None:
            break
        levels.append(g)
        if levels[-1].level != len(levels) - 1:  # a missing level number ends the hierarchy
            levels.pop()
            break
    if len(levels) == len(h0.levels) and all(g is g0 for g, g0 in zip(levels, h0.levels)):
        return h0 if extent == h0.extent else h0.with_extent(extent)
    return HierGeometry(levels, extent)


def locate(src, dst):
    """Store slot of the ``src`` patch holding every ``dst`` patch, or None if same layout.

    src, dst : LevelGeometry of the same level. Returns int64 ``(dst.npatch,)``. Raises
    ValueError when a dst patch is not inside a single src patch.
    """
    if src is dst or src.key() == dst.key():
        return None
    pairs = src.bin_index().pairs(dst.lower, dst.upper)
    inside = (pairs.lo == dst.lower[pairs.a]).all(axis=1) & (
        pairs.hi == dst.upper[pairs.a]).all(axis=1)
    slot = np.full(dst.npatch, -1, dtype=np.int64)
    slot[pairs.a[inside]] = pairs.b[inside]
    if (slot < 0).any():
        k = int(np.nonzero(slot < 0)[0][0])
        raise ValueError(
            f"level {dst.level}: target patch with cells {dst.lower[k]}..{dst.upper[k]} is not "
            "inside one patch of the field; restrict onto a layout made from this field "
            "(e.g. the result of an operation with it)"
        )
    return slot


def slot_blocks(sizes, max_nodes=None):
    """Consecutive slot ranges ``(i0, i1)`` holding at most ``max_nodes`` nodes each.

    max_nodes defaults to the module's ``BLOCK_NODES`` (read at call time). A single patch
    larger than ``max_nodes`` forms a block alone. One iteration per block.
    """
    max_nodes = BLOCK_NODES if max_nodes is None else max_nodes
    sizes = np.asarray(sizes, dtype=np.int64)
    ends = np.cumsum(sizes)
    n = len(sizes)
    i0 = 0
    while i0 < n:
        before = ends[i0] - sizes[i0]
        i1 = max(int(np.searchsorted(ends, before + max_nodes, side="right")), i0 + 1)
        yield i0, i1
        i0 = i1


def gather_indices(src, dst, slot, i0, i1):
    """Flat indices in ``src`` buffers of the nodes of dst patches ``i0 .. i1-1``.

    Order = dst flat order of those patches (x fastest per patch), so that
    ``dst_flat[dst.offset[i0]:dst.offset[i1]] = src_flat[idx]``.
    """
    s = slot[i0:i1]
    return box_flat_indices(src.offset[s], src.shape[s], dst.lower[i0:i1] - src.lower[s],
                            dst.shape[i0:i1])


def relayout_level(ld, dst):
    """LevelData of ``ld`` on the layout ``dst`` (shares buffers when the layout is the same)."""
    slot = locate(ld.geom, dst)
    if slot is None:
        return ld if ld.geom is dst else LevelData.from_arrays(dst, ld.time, ld.buffers())
    bufs = ld.buffers()
    out = {n: np.empty(dst.nnodes, dtype=ld.dtype) for n in ld.names}
    for i0, i1 in slot_blocks(dst.size):
        idx = gather_indices(ld.geom, dst, slot, i0, i1)
        n0, n1 = dst.offset[i0], dst.offset[i1]
        for n in ld.names:
            np.take(bufs[n], idx, out=out[n][n0:n1])
    return LevelData.from_arrays(dst, ld.time, out)


def _plan_check(snaps, hier, ncomp, dtypes, what):
    plan = ReadPlan()
    for s in snaps:
        for ld in s.levels[: len(hier.levels)]:
            plan.add_read(ld)
    nnodes = sum(g.nnodes for g in hier.levels)
    plan.add_result(nnodes, ncomp, *dtypes, what=what)
    plan.check(what=what)


def _result_dtype_args(operands):
    out = []
    for o in operands:
        if isinstance(o, Snapshot):
            out.append(o.levels[0].dtype)
        elif o is not None:
            out.append(o)  # Python scalars stay 'weak' (NEP 50)
    return out


# --------------------------------------------------------------------------- generic apply


# DESIGN-NOTE: operands at different times give a result stamped with the first field
# operand's time (a LevelData carries one time).
def apply(func, operands, what="operation"):
    """``func(*values)`` elementwise over operands (Snapshots or scalars) -> Snapshot.

    func : a numpy ufunc or any callable taking 1D arrays/scalars and returning a new 1D
           array of the same length (elementwise).
    Same layout: one ``func`` call per level and per component on whole flat buffers.
    Different layouts: per block of result patches, every operand is gathered, then one
    call. The result's time is the first field operand's time.
    """
    operands = [_check_operand(o, what) if o is not None else None for o in operands]
    snaps = _snapshots(operands)
    if not snaps:
        raise TypeError(f"{what}: at least one operand must be an AMR field")
    names = result_names(snaps, what)
    hier = common_hier([s.hier for s in snaps])
    _plan_check(snaps, hier, len(names), _result_dtype_args(operands), what)
    time = snaps[0].time
    levels = []
    for k, geom in enumerate(hier.levels):
        sources = []  # per operand: (kind, payload)
        for o in operands:
            if isinstance(o, Snapshot):
                ld = o.levels[k]
                sources.append(("field", ld, locate(ld.geom, geom)))
            else:
                sources.append(("scalar", o, None))
        for s in sources:  # every component an operand needs, in one read per level
            if s[0] == "field":
                s[1].ensure(None, sorted({_component_of(s[1], n) for n in names}))
        if all(src[2] is None for src in sources):
            out = {}
            for name in names:
                args = [s[1].buffer(_component_of(s[1], name)) if s[0] == "field" else s[1]
                        for s in sources]
                out[name] = np.asarray(func(*args))
        else:
            out = _apply_blocks(func, sources, geom, names)
        levels.append(LevelData.from_arrays(geom, time, out))
    return Snapshot(hier, tuple(levels))


def _apply_blocks(func, sources, geom, names):
    """Different layouts: gather every operand block by block, then call ``func``."""
    out = {}
    bufs = [
        {n: s[1].buffer(_component_of(s[1], n)) for n in names} if s[0] == "field" else None
        for s in sources
    ]
    for i0, i1 in slot_blocks(geom.size):
        n0, n1 = int(geom.offset[i0]), int(geom.offset[i1])
        idx = [gather_indices(s[1].geom, geom, s[2], i0, i1)
               if s[0] == "field" and s[2] is not None else None for s in sources]
        for name in names:
            args = []
            for s, b, ix in zip(sources, bufs, idx):
                if s[0] == "scalar":
                    args.append(s[1])
                elif ix is None:
                    args.append(b[name][n0:n1])
                else:
                    args.append(b[name][ix])
            r = np.asarray(func(*args))
            if name not in out:
                out[name] = np.empty(geom.nnodes, dtype=r.dtype)
            out[name][n0:n1] = r
    return out


# --------------------------------------------------------------------------- public API


def binary(ufunc, a, b):
    """``ufunc(a, b)`` elementwise; a, b: Snapshot or scalar (at least one Snapshot).

    ufunc : numpy ufunc with two inputs (np.add, np.multiply, np.power, np.greater, ...),
            or any elementwise callable ``f(x, y)`` on 1D arrays.
    Returns a Snapshot of in-memory data on the common layout (module doc). Reads the
    operands' missing patches first (memory check). Raises TypeError for arrays with
    ndim > 0, unknown operand types, or incompatible component names; ValueError for
    different grids or non-overlapping extents.
    """
    return apply(ufunc, [a, b], what=f"{getattr(ufunc, '__name__', 'binary')}")


def unary(ufunc, a, **kw):
    """``ufunc(a, **kw)`` elementwise on every component of Snapshot ``a`` (np.sqrt, np.abs ...).

    ``out=`` and ``where=`` are refused (TypeError): results are always new arrays.
    One ufunc call per level and per component on the whole flat buffers.
    """
    if not isinstance(a, Snapshot):
        raise TypeError(f"unary: expected an AMR field, got {type(a).__name__}")
    bad = [k for k in ("out", "where") if k in kw]
    if bad:
        raise TypeError(f"{bad[0]}= is not supported on AMR fields; use the returned field")
    what = getattr(ufunc, "__name__", "unary")
    _plan_check([a], a.hier, len(a.names), [a.levels[0].dtype], what)
    levels = tuple(
        LevelData.from_arrays(ld.geom, ld.time,
                              {n: np.asarray(ufunc(b, **kw)) for n, b in ld.buffers().items()})
        for ld in a.levels
    )
    return Snapshot(a.hier, levels)


def align(a, b):
    """Bring two Snapshots onto their common layout.

    Returns ``(hier, lds_a, lds_b)``: the common HierGeometry and, for each operand, a tuple
    of LevelData on it (the operand's own LevelData where the level is unchanged, gathered
    in-memory copies elsewhere). Raises ValueError for different grids.
    """
    for x in (a, b):
        if not isinstance(x, Snapshot):
            raise TypeError(f"align: expected two AMR fields, got {type(x).__name__}")
    hier = common_hier([a.hier, b.hier])
    la = tuple(relayout_level(ld, g) for ld, g in zip(a.levels, hier.levels))
    lb = tuple(relayout_level(ld, g) for ld, g in zip(b.levels, hier.levels))
    return hier, la, lb


def restrict(a, hier):
    """Snapshot ``a`` re-laid onto ``hier`` (whose extent the result takes).

    Every patch of ``hier`` at level l must lie inside one patch of ``a`` at level l (true
    for any layout produced by ``common_hier`` / ``binary`` with ``a``). Levels of ``hier``
    must all exist in ``a``. Raises ValueError otherwise, or for different grids.
    """
    if not isinstance(a, Snapshot):
        raise TypeError(f"restrict: expected an AMR field, got {type(a).__name__}")
    check_same_grid([a.hier, hier])
    levels = []
    for g in hier.levels:
        if not a.hier.has_level(g.level):
            raise ValueError(
                f"restrict: the target has level {g.level}, the field does not; restrict "
                "onto a layout made from this field"
            )
        k = a.hier.level_numbers().index(g.level)
        levels.append(relayout_level(a.levels[k], g))
    return Snapshot(hier, tuple(levels))


def where(cond, a, b):
    """``np.where(cond, a, b)`` elementwise; each argument a Snapshot or a scalar.

    ``cond`` is usually a boolean scalar field (e.g. ``binary(np.greater, n, 1.0)``); it is
    broadcast over vector components. Same layout and dtype rules as ``binary``.
    """
    return apply(np.where, [cond, a, b], what="where")


def clip(a, lo, hi):
    """``np.clip(a, lo, hi)`` elementwise; lo/hi are scalars, fields, or None (no bound)."""
    if lo is None and hi is None:
        raise ValueError("clip: give at least one of lo, hi")
    if lo is None:
        return apply(lambda x, h: np.clip(x, None, h), [a, hi], what="clip")
    if hi is None:
        return apply(lambda x, low: np.clip(x, low, None), [a, lo], what="clip")
    return apply(np.clip, [a, lo, hi], what="clip")


def component(a, name):
    """Component ``name`` of Snapshot ``a`` as a scalar Snapshot ``("value",)``, no copy.

    The flat buffers are shared (file-backed data is read first, then shared read-only).
    Raises KeyError naming the available components.
    """
    if not isinstance(a, Snapshot):
        raise TypeError(f"component: expected an AMR field, got {type(a).__name__}")
    if name not in a.names:
        raise KeyError(f"no component {name!r}; this field has {a.names}")
    levels = tuple(LevelData.from_arrays(ld.geom, ld.time, {"value": ld.buffer(name)})
                   for ld in a.levels)
    return Snapshot(a.hier, levels)


def stack_components(snaps, names):
    """One Snapshot whose components ``names`` are the scalar Snapshots ``snaps``.

    Example: ``stack_components([bx, by, bz], ("x", "y", "z"))``. Buffers are shared when
    every input is on the same layout; otherwise each one is gathered onto the common
    layout. Mixed dtypes are cast to ``np.result_type`` (a copy for the cast ones).
    Raises TypeError if an input is not a scalar field or the counts differ.
    """
    snaps, names = list(snaps), tuple(names)
    if len(snaps) != len(names) or not snaps:
        raise TypeError(f"stack_components: {len(snaps)} fields for names {names}")
    for s in snaps:
        if not isinstance(s, Snapshot) or s.names != SCALAR_NAMES:
            raise TypeError("stack_components: every input must be a scalar field ('value',)")
    hier = common_hier([s.hier for s in snaps])
    dtype = np.result_type(*[s.levels[0].dtype for s in snaps])
    levels = []
    for k, g in enumerate(hier.levels):
        lds = [relayout_level(s.levels[k], g) for s in snaps]
        arrays = {n: ld.buffer("value").astype(dtype, copy=False) for n, ld in zip(names, lds)}
        levels.append(LevelData.from_arrays(g, snaps[0].time, arrays))
    return Snapshot(hier, tuple(levels))


# --------------------------------------------------------------------------- vector algebra


def _require_vector(x, what):
    if not isinstance(x, Snapshot) or tuple(x.names) != VECTOR_NAMES:
        got = x.names if isinstance(x, Snapshot) else type(x).__name__
        raise TypeError(f"{what} expects vector fields with components {VECTOR_NAMES}, got {got}")


def _blockwise(func, arrays, dtype):
    """``func(*[a[block] for a in arrays])`` written block by block into a new array.

    Bounds the temporaries of a compound expression to BLOCK_NODES values each.
    """
    n, step = len(arrays[0]), BLOCK_NODES
    out = np.empty(n, dtype=dtype)
    for c in range(0, n, step):
        sl = slice(c, min(c + step, n))
        out[sl] = func(*[a[sl] for a in arrays])
    return out


def dot(a, b):
    """Scalar product ``a.x b.x + a.y b.y + a.z b.z`` of two vector Snapshots."""
    _require_vector(a, "dot")
    _require_vector(b, "dot")
    hier, la, lb = align(a, b)
    dtype = np.result_type(a.levels[0].dtype, b.levels[0].dtype)
    _plan_check([], hier, 1, [dtype], "dot")
    levels = []
    for g, x, y in zip(hier.levels, la, lb):
        X = list(x.buffers(VECTOR_NAMES).values())
        Y = list(y.buffers(VECTOR_NAMES).values())
        r = _blockwise(lambda x0, x1, x2, y0, y1, y2: x0 * y0 + x1 * y1 + x2 * y2,
                       X + Y, dtype)
        levels.append(LevelData.from_arrays(g, a.time, {"value": r}))
    return Snapshot(hier, tuple(levels))


def cross(a, b):
    """Vector product ``a x b`` of two vector Snapshots -> vector Snapshot."""
    _require_vector(a, "cross")
    _require_vector(b, "cross")
    hier, la, lb = align(a, b)
    dtype = np.result_type(a.levels[0].dtype, b.levels[0].dtype)
    _plan_check([], hier, 3, [dtype], "cross")
    levels = []
    for g, x, y in zip(hier.levels, la, lb):
        X = list(x.buffers(VECTOR_NAMES).values())
        Y = list(y.buffers(VECTOR_NAMES).values())
        out = {}
        for i, c in enumerate(VECTOR_NAMES):
            j, k = (i + 1) % 3, (i + 2) % 3
            out[c] = _blockwise(lambda xj, xk, yj, yk: xj * yk - xk * yj,
                                [X[j], X[k], Y[j], Y[k]], dtype)
        levels.append(LevelData.from_arrays(g, a.time, out))
    return Snapshot(hier, tuple(levels))


def modulus(a):
    """``sqrt(a.x**2 + a.y**2 + a.z**2)`` of a vector Snapshot -> scalar Snapshot.

    dtype: float of the input (float32 stays float32).
    """
    _require_vector(a, "modulus")
    dtype = np.result_type(a.levels[0].dtype, np.float16)
    _plan_check([a], a.hier, 1, [dtype], "modulus")
    levels = []
    for ld in a.levels:
        X = list(ld.buffers(VECTOR_NAMES).values())
        r = _blockwise(lambda x, y, z: np.sqrt(x * x + y * y + z * z), X, dtype)
        levels.append(LevelData.from_arrays(ld.geom, ld.time, {"value": r}))
    return Snapshot(a.hier, tuple(levels))
