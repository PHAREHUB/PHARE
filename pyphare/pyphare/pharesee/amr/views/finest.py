"""``finest()``: an AMR snapshot painted onto one regular node grid.

The output is the node grid of one level ``L`` (the finest level asked for) over a box.
It is built level by level, coarse to fine:

1. level 0: the patches that touch the level-0 *window* (the box coarsened to level 0,
   plus one node of margin, clipped to the domain) are copied into that window;
2. the window is upsampled by 2 into the window of the next level: nodes that exist on both
   grids are copied, the new nodes are filled (``fill="cubic"``: clipped Catmull-Rom along each
   axis in turn, linear next to the window edges; ``fill="nearest"``: the coarse node below);
3. that level's patches are copied on top, and so on; the last upsample and the last
   paint write straight into the output.

So every output node that is a node of some level holding it in a patch shows that
level's value, bit for bit (the finest such level wins); the other nodes ("gaps") are
display interpolation of coarser data, never a new extremum. A NaN in the coarse data
under a fine patch (e.g. a masked divB node) is hidden where the fine level paints, but can
reach the gap nodes next to it through the interpolation stencil. ``level_map`` says which level
each node shows (for a gap: the coarser of the two coarse nodes it lies between).

Memory: the output, its ``level_map`` (1 byte per node), the window of level ``L - 1`` (a
quarter of the output in 2D) and bounded slabs of scratch. Upsampling goes slab by slab
along the last axis, so no full-size temporary is ever made.
"""

from __future__ import annotations

import math
import time as _time

import numpy as np

from .. import plan
from ..boxes import AXES, PhysBox

#: level_map value of nodes no level-0 patch holds (data there is NaN)
NO_DATA = 255
#: target bytes of one output slab during upsampling (scratch is a few times this)
SLAB_BYTES = 8 * 2**20
FILLS = ("cubic", "nearest")


# --------------------------------------------------------------------------- result


class FinestGrid:
    """A snapshot painted on one regular node grid (see ``finest``).

    Attributes
      data      ndarray indexed ``[x, y(, z)]`` (x fastest in memory) for a one-component
                snapshot; dict ``name -> ndarray`` otherwise
      x, y, z   1D float64 node coordinates along each axis (``y`` absent in 1D, ``z`` only
                in 3D); node ``i`` is at ``(first_node + i) * dl``
      extent    pixel-edge extent for ``imshow``: ``(x[0] - dx/2, x[-1] + dx/2, y[0] - dy/2,
                y[-1] + dy/2)`` in 2D (one pair per axis in 1D/3D)
      box       PhysBox spanned by the output nodes (the asked box snapped outward to nodes)
      dl        (dim,) float64 spacing of the output grid
      time      time of the snapshot
      level     the level whose node grid this is
      level_map uint8 array, same shape as the data: the AMR level each node shows
                (``NO_DATA`` = 255 where no patch holds data; the value there is NaN)
    Unpacks as ``data, coords = grid`` with ``coords = (x, y[, z])``.
    """

    def __init__(self, data, coords, dl, time, level, level_map, single):
        self.data = data
        self._coords = tuple(coords)
        self.dl = np.asarray(dl, dtype=np.float64)
        self.time = time
        self.level = level
        self.level_map = level_map
        self._single = single
        self._timings = {}

    @property
    def dim(self):
        return len(self._coords)

    @property
    def x(self):
        return self._coords[0]

    @property
    def y(self):
        if self.dim < 2:
            raise AttributeError("a 1D grid has no y coordinate")
        return self._coords[1]

    @property
    def z(self):
        if self.dim < 3:
            raise AttributeError(f"a {self.dim}D grid has no z coordinate")
        return self._coords[2]

    @property
    def coords(self):
        """Tuple of the node coordinate arrays ``(x, y[, z])``."""
        return self._coords

    @property
    def names(self):
        """Component names held."""
        return (self._single,) if self._single is not None else tuple(self.data)

    @property
    def shape(self):
        return self.level_map.shape

    @property
    def box(self):
        return PhysBox(tuple(float(c[0]) for c in self._coords),
                       tuple(float(c[-1]) for c in self._coords))

    # DESIGN-NOTE: extent is the PIXEL-EDGE extent (nodes are pixel centres, half a cell
    # beyond the first/last node) so that imshow(data.T, extent=g.extent) puts every node at
    # its coordinate; the node bounds are g.box.
    @property
    def extent(self):
        out = []
        for c, h in zip(self._coords, self.dl):
            out += [float(c[0] - h / 2), float(c[-1] + h / 2)]
        return tuple(out)

    def __getitem__(self, name):
        """Array of component ``name`` (for a one-component grid, its own name)."""
        if self._single is not None:
            if name == self._single:
                return self.data
            raise KeyError(f"this grid holds one component, {self._single!r}; use g.data")
        try:
            return self.data[name]
        except KeyError:
            raise KeyError(f"unknown component {name!r}; this grid has {self.names}") from None

    def __iter__(self):
        """``data, coords = grid``."""
        return iter((self.data, self._coords))

    def __array__(self, dtype=None, copy=None):
        if self._single is None:
            raise TypeError(
                f"this grid has {len(self.data)} components {self.names}: take one with "
                f"g['{self.names[0]}'] before converting to an array"
            )
        a = self.data
        if dtype is not None and np.dtype(dtype) != a.dtype:
            return a.astype(dtype)
        return a.copy() if copy else a

    def imshow(self, ax=None, name=None, **kw):
        """Draw a 2D grid with matplotlib ``imshow`` and return the ``AxesImage``.

        Transposes the ``[x, y]`` data, sets ``origin="lower"``, the pixel-edge ``extent``
        and ``interpolation="nearest"`` (each overridable through ``kw``).
        ax   : matplotlib Axes (default: ``plt.gca()``)
        name : component to draw (required for a multi-component grid)
        Raises ValueError for 1D/3D grids (cut or reduce them first).
        """
        if self.dim != 2:
            raise ValueError(f"imshow needs a 2D grid, this one is {self.dim}D; "
                             "take a 2D slice of g.data and draw it yourself")
        if name is None:
            if self._single is None:
                raise ValueError(f"name= is needed: this grid has components {self.names}")
            a = self.data
        else:
            a = self[name]
        if ax is None:
            import matplotlib.pyplot as plt

            ax = plt.gca()
        opts = dict(origin="lower", extent=self.extent, interpolation="nearest")
        opts.update(kw)
        return ax.imshow(a.T, **opts)

    def __repr__(self):
        shape = "x".join(str(n) for n in self.shape)
        return (f"FinestGrid({self.dim}D {shape} nodes at level {self.level}, t={self.time:g}, "
                f"components={self.names}, box={self.box})")


# --------------------------------------------------------------------------- upsampling


def _up_axis(src, dst, axis, f0, cubic):
    """Ratio-2 upsampling of ``src`` along ``axis`` into ``dst``.

    ``dst`` node ``j`` along ``axis`` is fine node ``k = f0 + j`` where coarse node ``i`` of
    ``src`` is fine node ``2 i``. Even ``k``: copy of coarse ``k / 2``. Odd ``k``: between
    coarse ``b = (k-1)/2`` and ``c = b + 1``: cubic -> Catmull-Rom
    ``(9 (f_b + f_c) - f_{b-1} - f_{c+1}) / 16`` clipped to ``[min, max](f_b, f_c)``, or
    ``(f_b + f_c) / 2`` when ``b - 1`` or ``c + 1`` is outside ``src``; nearest -> ``f_b``.
    Other axes of ``src`` and ``dst`` must have equal sizes. Vectorised slices, no loops.
    """
    s = np.moveaxis(src, axis, 0)
    d = np.moveaxis(dst, axis, 0)
    nc, n = s.shape[0], d.shape[0]
    je = f0 % 2  # first dst position holding an even fine node
    ce = (f0 + je) // 2
    ne = (n - je + 1) // 2
    if ne > 0:
        d[je::2] = s[ce : ce + ne]
    jo = 1 - je
    no = (n - jo + 1) // 2
    if no <= 0:
        return
    co = (f0 + jo - 1) // 2  # coarse node left of the first odd fine node
    b = s[co : co + no]
    c = s[co + 1 : co + 1 + no]
    m = d[jo::2]
    if not cubic:
        m[...] = b
        return
    # Catmull-Rom needs b-1 and c+1: interior odd nodes are i in [i0, i1)
    i0 = min(max(0, 1 - co), no)
    i1 = max(min(no, nc - 2 - co), i0)
    half = m.dtype.type(0.5)
    for lo, hi in ((0, i0), (i1, no)):  # at most one node each, next to the window edges
        if hi > lo:
            np.add(b[lo:hi], c[lo:hi], out=m[lo:hi])
            m[lo:hi] *= half
    if i1 <= i0:
        return
    bi, ci, mi = b[i0:i1], c[i0:i1], m[i0:i1]
    np.add(bi, ci, out=mi)
    mi *= m.dtype.type(9)
    mi -= s[co + i0 - 1 : co + i1 - 1]
    mi -= s[co + i0 + 2 : co + i1 + 2]
    mi *= m.dtype.type(1 / 16)
    # no new extrema: clip to the two neighbouring coarse values
    np.maximum(mi, np.minimum(bi, ci), out=mi)
    np.minimum(mi, np.maximum(bi, ci), out=mi)


# DESIGN-NOTE: DESIGN does not say which level a gap node 'shows'; we take the coarser of the
# two coarse nodes it lies between (cubic) or the node it copies (nearest). Nodes no level-0
# patch holds get NO_DATA (255) and NaN.
def _up_axis_map(src, dst, axis, f0, cubic):
    """``_up_axis`` for the level map: even nodes copied, odd nodes get the coarser
    (smaller) level of their two neighbours (``fill="cubic"``) or the lower one (nearest)."""
    s = np.moveaxis(src, axis, 0)
    d = np.moveaxis(dst, axis, 0)
    n = d.shape[0]
    je = f0 % 2
    ce = (f0 + je) // 2
    ne = (n - je + 1) // 2
    if ne > 0:
        d[je::2] = s[ce : ce + ne]
    jo = 1 - je
    no = (n - jo + 1) // 2
    if no <= 0:
        return
    co = (f0 + jo - 1) // 2
    if cubic:
        np.minimum(s[co : co + no], s[co + 1 : co + 1 + no], out=d[jo::2])
    else:
        d[jo::2] = s[co : co + no]


def upsample(src, dst, f0, cubic, kernel=_up_axis, slab_bytes=None):
    """Write the ratio-2 upsampling of window ``src`` into window ``dst`` (both ``[x, y, z]``).

    f0 : (dim,) int, fine index of ``dst[0, 0, 0]`` counted from ``2 * (coarse index of
         src[0, 0, 0])``. ``src`` must cover every coarse node the fill needs (one node of
         margin beyond the fine range, except at the window/domain edges).
    Works slab by slab along the last axis: for each slab of ``dst``, the coarse rows it
    needs (+1 margin) are upsampled along axis 0, then 1, ..., the last pass writing into
    ``dst`` directly. Scratch per slab is about ``slab_bytes`` (default ``SLAB_BYTES``).
    """
    dim = dst.ndim
    last = dim - 1
    f0 = [int(v) for v in f0]
    slab_bytes = SLAB_BYTES if slab_bytes is None else slab_bytes
    row = math.prod(dst.shape[:last]) * dst.itemsize
    rows = max(4, int(slab_bytes // max(row, 1)))
    nlast, nc = dst.shape[last], src.shape[last]
    for s0 in range(0, nlast, rows):
        s1 = min(s0 + rows, nlast)
        k0, k1 = f0[last] + s0, f0[last] + s1 - 1  # fine range along the last axis
        c0 = max(k0 // 2 - 1, 0)
        c1 = min(-(-k1 // 2) + 1, nc - 1)
        cur = src[..., c0 : c1 + 1]
        for ax in range(last):  # axes before the last: into bounded scratch
            shape = list(cur.shape)
            shape[ax] = dst.shape[ax]
            tmp = np.empty(shape, dtype=dst.dtype, order="F")
            kernel(cur, tmp, ax, f0[ax], cubic)
            cur = tmp
        kernel(cur, dst[..., s0:s1], last, k0 - 2 * c0, cubic)


# --------------------------------------------------------------------------- painting


def _flat(ld, name):
    """Flat buffer of component ``name`` of a LevelData whose needed patches were ensured.

    ``LevelData.raw``: never reads (``buffer`` would read every missing patch first); valid
    at the ensured patches only.
    """
    return ld.raw(name)


class _LevelPaint:
    """Patches of one level that touch one window, with their copy slices (plain ints)."""

    def __init__(self, geom, wlo, whi):
        self.geom = geom
        dim = geom.dim
        # node box [lower, upper + 1] meets [wlo, whi]  <=>  cells meet [wlo - 1, whi]
        clo = np.clip(np.asarray(wlo) - 1, 0, geom.domain_cells - 1)
        chi = np.clip(np.asarray(whi), 0, geom.domain_cells - 1)
        slots = geom.patches_touching(clo, chi)
        lower = geom.lower[slots]
        plo = np.maximum(lower, wlo)
        phi = np.minimum(geom.upper[slots] + 1, whi)
        keep = (plo <= phi).all(axis=1)
        self.slots = slots[keep]
        plo, phi, lower = plo[keep], phi[keep], lower[keep]
        # rows of (window start, patch start, length) per axis, as Python ints
        self.rows = np.concatenate([plo - wlo, plo - lower, phi - plo + 1], axis=1).tolist()
        self.dim = dim

    def paint(self, flat, win, value_map=None, level=0):
        """Copy every patch's part inside the window (one slice copy per patch)."""
        geom, dim = self.geom, self.dim
        offset, shape = geom.offset, geom.shape
        for slot, r in zip(self.slots.tolist(), self.rows):
            w = tuple(slice(r[a], r[a] + r[2 * dim + a]) for a in range(dim))
            if flat is not None:
                p = tuple(slice(r[dim + a], r[dim + a] + r[2 * dim + a]) for a in range(dim))
                v = flat[offset[slot] : offset[slot + 1]].reshape(tuple(shape[slot]), order="F")
                win[w] = v[p]
            if value_map is not None:
                value_map[w] = level


# --------------------------------------------------------------------------- planning


def _snap(q):
    r = np.round(q)
    return np.where(np.abs(q - r) <= 1e-9 * np.maximum(1.0, np.abs(q)), r, q)


def _output_nodes(box, dlL, NL):
    """Node range ``[lo, hi]`` (int64 (dim,)) of level-L nodes covering ``box`` (snapped out)."""
    lo, hi = box.clipped(NL * dlL)
    if (lo > hi).any():
        a = int(np.nonzero(lo > hi)[0][0])
        raise ValueError(f"box {box.to_dict()} lies outside the domain [0, {NL[a] * dlL[a]}] "
                         f"along {AXES[a]}; choose bounds inside the domain")
    nlo = np.floor(_snap(lo / dlL)).astype(np.int64)
    nhi = np.ceil(_snap(hi / dlL)).astype(np.int64)
    return np.clip(nlo, 0, NL), np.clip(nhi, 0, NL)


# DESIGN-NOTE: windows are clipped to the domain, never wrapped on periodic axes: gap nodes in
# the last coarse cell next to a periodic boundary are filled linearly, not by Catmull-Rom.
def _windows(nlo, nhi, N0, L):
    """Node windows ``[(lo, hi)]`` for levels 0..L: level L = the output, coarser ones the
    finer window coarsened plus one node of margin, clipped to the domain."""
    wins = [(nlo, nhi)]
    for l in range(L - 1, -1, -1):
        flo, fhi = wins[0]
        Nl = N0 * 2**l
        lo = np.maximum(flo // 2 - 1, 0)
        hi = np.minimum(-(-fhi // 2) + 1, Nl)
        wins.insert(0, (lo, hi))
    return wins


# DESIGN-NOTE: a box given together with an extent is intersected with it (empty -> ValueError).
def _resolve_box(hier, box):
    box = PhysBox.parse(box, hier.dim)
    if hier.extent is not None:
        if box is None:
            return hier.extent
        inside = box.intersect(hier.extent)
        if inside is None:
            raise ValueError(f"box {box} does not meet the extent of this field "
                             f"{hier.extent}; choose a box inside it")
        return inside
    return box if box is not None else PhysBox((0.0,) * hier.dim, tuple(hier.domain))


def _names(snap, names):
    all_names = tuple(snap.names)
    if names is None:
        return all_names
    wanted = (names,) if isinstance(names, str) else tuple(names)
    unknown = [n for n in wanted if n not in all_names]
    if unknown or not wanted:
        raise KeyError(f"unknown component(s) {unknown or wanted}; this field has {all_names}")
    return wanted


# --------------------------------------------------------------------------- finest


def finest(snap, box=None, max_level=None, fill="cubic", dtype=None, force=False, names=None):
    """Paint a snapshot on the regular node grid of its finest level over a box.

    snap      : ``store.Snapshot`` (one quantity at one time)
    box       : ``dict(x=(x0, x1), y=(y0, y1))`` / ``PhysBox`` / None, physical units. None ->
                the snapshot's extent, or the whole domain. A box is intersected with the
                extent; it is snapped OUTWARD to the output level's nodes and clipped to the
                domain (a box narrower than a cell gives the two nodes around it).
    max_level : resolution level of the output: ``L = min(max_level, finest level present)``.
                Levels above L are never read.
    fill      : how nodes that no level holds are filled: "cubic" (default; Catmull-Rom along
                each axis in turn, clipped to the two neighbouring coarse values, linear at
                the window edges) or "nearest" (coarse node below, ties toward the lower
                index). These values are display interpolation, not data.
    dtype     : output dtype (default: the field's); must be a float type.
    force     : go ahead when the memory check says it does not fit.
    names     : component name or list of names to paint (default: all).

    Returns a ``FinestGrid``: ``.data[x, y(, z)]`` (or dict per component), node coordinates
    ``.x .y .z``, ``.level_map``. Fine values are bit-identical to the patch data (when
    ``dtype`` is the field's). Reads exactly the patches of levels 0..L that touch the
    windows (box coarsened per level, one node of margin): one ``ensure`` call per level.
    Raises PhareseeMemoryError when output + scratch do not fit (unless ``force``),
    ValueError for a box outside the domain/extent, a bad ``fill`` or ``max_level``,
    KeyError for unknown ``names``, TypeError for a non-float dtype.
    """
    t_start = _time.perf_counter()
    if fill not in FILLS:
        raise ValueError(f"fill={fill!r}: choose one of {FILLS}")
    cubic = fill == "cubic"
    hier = snap.hier
    lds = {g.level: ld for g, ld in zip(hier.levels, snap.levels)}
    finest_present = hier.finest_level
    if max_level is None:
        L = finest_present
    else:
        if int(max_level) != max_level or max_level < 0:
            raise ValueError(f"max_level={max_level!r} must be an integer >= 0")
        L = min(int(max_level), finest_present)
    sel = _names(snap, names)
    single = snap.names[0] if len(snap.names) == 1 else None
    out_dtype = np.dtype(dtype) if dtype is not None else np.dtype(snap.levels[0].dtype)
    # DESIGN-NOTE: dtype= must be a float type (gaps are interpolated, missing data is NaN).
    if out_dtype.kind != "f":
        raise TypeError(f"dtype={out_dtype} is not a float type; finest() fills gaps by "
                        "interpolation and marks missing data with NaN")

    g0 = hier.levels[0]
    N0 = np.asarray(g0.domain_cells, np.int64)
    dl0 = np.asarray(g0.dl, np.float64)
    pbox = _resolve_box(hier, box)
    dlL = dl0 / 2**L
    nlo, nhi = _output_nodes(pbox, dlL, N0 * 2**L)
    wins = _windows(nlo, nhi, N0, L)
    out_shape = tuple(int(v) for v in nhi - nlo + 1)
    nout = math.prod(out_shape)

    # patches to paint, per level (geometry only, no data yet)
    paints = {}
    for l in range(L + 1):
        if l in lds:
            paints[l] = _LevelPaint(lds[l].geom, *wins[l])

    # memory: outputs + level map + the largest pair of consecutive windows + slabs + reads
    win_nodes = [math.prod(int(v) for v in hi - lo + 1) for lo, hi in wins]
    scratch = 0
    for l in range(L):
        scratch = max(scratch, win_nodes[l] * (out_dtype.itemsize + 1)
                      + (win_nodes[l + 1] * (out_dtype.itemsize + 1) if l + 1 < L else 0))
    need = plan.ReadPlan()
    need.add(nout * out_dtype.itemsize * len(sel), f"finest() output {out_shape}")
    need.add(nout, "level map")
    need.add(scratch + 4 * SLAB_BYTES, "coarser windows and slabs")
    for l, p in paints.items():
        need.add_read(lds[l], p.slots, sel)
    need.check(force=force, what=f"finest() on a {'x'.join(map(str, out_shape))} grid")

    t_read = _time.perf_counter()
    for l, p in paints.items():
        if len(p.slots):
            lds[l].ensure(p.slots, sel, force=force)
    timings = dict(plan=t_read - t_start, read=_time.perf_counter() - t_read, paint=0.0,
                   upsample=0.0, map=0.0)

    # DESIGN-NOTE: components share the plan (windows, patch lists, slices) and the level map;
    # the coarse windows are re-allocated per component and freed level by level, which keeps
    # the peak at output + one window instead of holding every window across components.
    level_map = None
    missing0 = None  # level-0 window nodes no level-0 patch holds (None: there are none)
    data = {}
    for name in sel:
        need_map = level_map is None
        win, wmap = None, None
        for l in range(L + 1):
            shape = tuple(int(v) for v in wins[l][1] - wins[l][0] + 1)
            new = np.empty(shape, dtype=out_dtype, order="F")
            nmap = np.empty(shape, dtype=np.uint8, order="F") if need_map else None
            t = _time.perf_counter()
            if win is None:
                if nmap is not None:
                    nmap.fill(NO_DATA)
            else:
                f0 = wins[l][0] - 2 * wins[l - 1][0]
                upsample(win, new, f0, cubic)
                timings["upsample"] += _time.perf_counter() - t
                t = _time.perf_counter()
                if nmap is not None:
                    upsample(wmap, nmap, f0, cubic, kernel=_up_axis_map)
                timings["map"] += _time.perf_counter() - t
            win, wmap = None, None  # free the coarser window before painting
            t = _time.perf_counter()
            p = paints.get(l)
            if p is not None and len(p.slots):
                p.paint(_flat(lds[l], name), new, nmap, l)
            if l == 0:
                if need_map:
                    missing0 = nmap == NO_DATA
                    if not missing0.any():
                        missing0 = None
                if missing0 is not None:
                    new[missing0] = np.nan
            timings["paint"] += _time.perf_counter() - t
            win, wmap = new, nmap
        data[name] = win
        if need_map:
            level_map = wmap
    timings["total"] = _time.perf_counter() - t_start

    coords = [(nlo[a] + np.arange(out_shape[a], dtype=np.float64)) * dlL[a]
              for a in range(hier.dim)]
    if single is not None:
        grid = FinestGrid(data[single], coords, dlL, snap.time, L, level_map, single)
    else:
        grid = FinestGrid(data, coords, dlL, snap.time, L, level_map, None)
    grid._timings = timings
    return grid
