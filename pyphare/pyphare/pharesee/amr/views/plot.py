"""``field.plot(...)``: master's plotting keywords on the new fields (DESIGN 2 "plot()", 8.4).

2D: ONE ``imshow`` of ``field.finest(box)`` where the box comes from ``xlim``/``ylim`` (or the
field's extent), so a zoom reads only the patches under it. 1D: one line per level (patches
separated by gaps). 3D: ValueError asking for a cut.

Colour limits (2D), in this order of precedence:
  * explicit ``vmin``/``vmax`` (each one wins on its side) or ``norm=``;
  * ``robust=True``: percentiles 1 and 99.7 of the shown data;
  * default: min/max of the shown data (NaN ignored) - master's rule restricted to what is
    on screen (master took the global min/max over every level and ghost node, #1450);
  * ``symmetric=True``: limits made ``(-m, m)`` with ``m`` the largest absolute limit.
The colourbar gets ``extend`` arrows on every side where shown data lie beyond the limits.

matplotlib is imported inside the functions only.
"""

from __future__ import annotations

import numpy as np

from .finest import NO_DATA, FinestGrid

#: master's keywords (hierarchy.py plot1d 396-447 and plot2d 449-541)
MASTER_KEYWORDS = (
    "qty", "time", "levels", "ax", "cmap", "vmin", "vmax", "plot_patches", "patchcolors",
    "lw", "ls", "aspect", "title", "xlabel", "ylabel", "xlim", "ylim", "legend", "filename",
    "dpi", "marker", "color",
)
#: keywords added by the new plot
NEW_KEYWORDS = ("robust", "symmetric", "norm", "max_level", "fill")
#: keywords master's callers pass and master's plot never read (it took any **kwargs):
#: accepted and ignored (tests/simulator/test_run.py passes all_primal=False to J.plot)
IGNORED_MASTER_KEYWORDS = ("all_primal",)

ROBUST_PERCENTILES = (1.0, 99.7)
DEFAULT_CMAP = "Spectral_r"  # master's default (hierarchy.py:504)


def _check_keywords(kwargs):
    accepted = MASTER_KEYWORDS + NEW_KEYWORDS + IGNORED_MASTER_KEYWORDS
    unknown = sorted(k for k in kwargs if k not in accepted)
    if unknown:
        raise TypeError(
            f"plot() got unexpected keyword argument(s) {unknown}; accepted: "
            f"{', '.join(MASTER_KEYWORDS + NEW_KEYWORDS)}"
        )


def _one_time(field, time):
    if time is not None:
        return field.at(time)
    field.snapshot()  # ValueError listing the times when the field holds several
    return field


def _component(field, qty):
    """ScalarField to draw: the field itself (scalar) or its component ``qty``."""
    if tuple(field.names) == ("value",):
        # a component field answers its parent's name too (B.x.plot(qty="x"), as master)
        alias = getattr(field, "_component_alias", None)
        if qty in (None, "value") or qty == field.quantity or (alias is not None
                                                               and qty == alias):
            return field
        return field[qty]
    if qty is None:
        raise ValueError(
            "multiple quantities in the field, please specify a quantity in "
            + " ".join(field.names)
        )
    return field[qty]


def default_title(field):
    """Title used when ``title=`` is not given: "" for stored quantities; for a derived
    field with a stencil (J, divB), its name and stencil, e.g. ``"J_z (stencil: node)"``."""
    if getattr(field, "stencil", None) is None:
        return ""
    return f"{field.quantity or ''} (stencil: {field.stencil})".strip()


def _levels(levels, present, default_all):
    if levels is None:
        return tuple(present) if default_all else (0,)
    if np.ndim(levels) == 0:
        levels = (levels,)
    return tuple(sorted({int(v) for v in levels}))


def _per_level(option, levels, default):
    """master's per-level option: a dict, a sequence zipped with the levels, or one value."""
    if option is None:
        return dict.fromkeys(levels, default)
    if isinstance(option, dict):
        return {l: option.get(l, default) for l in levels}
    if isinstance(option, (list, tuple)):
        return dict(zip(levels, option))
    return dict.fromkeys(levels, option)


def _figure(ax):
    import matplotlib.pyplot as plt

    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.figure
    return fig, ax


def _limits(values, vmin, vmax, robust, symmetric):
    """Colour limits from the shown finite ``values`` (rules of the module docstring)."""
    lo = hi = None
    if values.size:
        if robust:
            lo, hi = (float(v) for v in np.percentile(values, ROBUST_PERCENTILES))
        else:
            lo, hi = float(values.min()), float(values.max())
        if symmetric:
            m = max(abs(lo), abs(hi))
            lo, hi = -m, m
    lo = lo if vmin is None else float(vmin)
    hi = hi if vmax is None else float(vmax)
    return lo, hi


# DESIGN-NOTE: the colourbar gets extend arrows whenever shown data lie beyond the limits,
# not only with robust=True (explicit vmin/vmax that clip the data are flagged too).
def _extend(values, lo, hi):
    """Colourbar ``extend`` for data ``values`` drawn with limits ``(lo, hi)``."""
    if not values.size or lo is None or hi is None:
        return "neither"
    below, above = bool(values.min() < lo), bool(values.max() > hi)
    return {(False, False): "neither", (True, False): "min", (False, True): "max",
            (True, True): "both"}[(below, above)]


def _view_box(kw):
    """PhysBox dict from xlim/ylim (None when neither is given)."""
    box = {}
    for axis, key in (("x", "xlim"), ("y", "ylim")):
        lim = kw.get(key)
        if lim is not None:
            a, b = float(lim[0]), float(lim[1])
            box[axis] = (min(a, b), max(a, b))
    return box or None


def patch_outline_segments(geom, view_lo, view_hi):
    """Outlines of the patches of one level that intersect a view, for a LineCollection.

    geom             : LevelGeometry (2D)
    view_lo, view_hi : (2,) physical bounds of the view
    Returns float64 ``(n, 5, 2)``: one closed polyline per patch whose physical extent
    ``[lower * dl, (upper + 1) * dl]`` overlaps the view with a non-zero area (a patch that
    only touches the view edge is left out). No per-patch Python loop.
    """
    lo = geom.lower * geom.dl
    hi = (geom.upper + 1) * geom.dl
    keep = ((hi > np.asarray(view_lo)) & (lo < np.asarray(view_hi))).all(axis=1)
    lo, hi = lo[keep], hi[keep]
    xs = np.stack([lo[:, 0], hi[:, 0], hi[:, 0], lo[:, 0], lo[:, 0]], axis=1)
    ys = np.stack([lo[:, 1], lo[:, 1], hi[:, 1], hi[:, 1], lo[:, 1]], axis=1)
    return np.stack([xs, ys], axis=-1)


def _plot2d(field, kw):
    from matplotlib.collections import LineCollection
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    snap = field.snapshot()
    present = [g.level for g in snap.hier.levels]
    levels = _levels(kw.get("levels"), present, default_all=True)
    max_level = max(levels) if kw.get("levels") is not None else None
    if kw.get("max_level") is not None:
        max_level = int(kw["max_level"]) if max_level is None else min(max_level,
                                                                       int(kw["max_level"]))
    g = field.finest(box=_view_box(kw), max_level=max_level, fill=kw.get("fill", "cubic"))
    data = np.asarray(g.data)
    shown = g.level_map != NO_DATA
    if kw.get("levels") is not None:
        shown &= np.isin(g.level_map, np.asarray(levels, dtype=np.uint8))
    finite = shown & np.isfinite(data)
    values = data[finite]

    norm = kw.get("norm")
    if norm is not None and (kw.get("vmin") is not None or kw.get("vmax") is not None):
        raise ValueError("plot(): give norm= or vmin=/vmax=, not both (set the limits on "
                         "the norm)")
    opts = dict(cmap=kw.get("cmap", DEFAULT_CMAP))
    if norm is not None:
        opts["norm"] = norm
    else:
        opts["vmin"], opts["vmax"] = _limits(values, kw.get("vmin"), kw.get("vmax"),
                                             kw.get("robust", False), kw.get("symmetric", False))

    fig, ax = _figure(kw.get("ax"))
    masked = np.ma.masked_array(data, mask=~finite)
    grid = FinestGrid(masked, g.coords, g.dl, g.time, g.level, g.level_map, g.names[0])
    im = grid.imshow(ax, **opts)
    extend = _extend(values, im.norm.vmin, im.norm.vmax)

    if kw.get("plot_patches", False) is True:
        colors = _per_level(kw.get("patchcolors"), levels, "k")
        widths = _per_level(kw.get("lw"), levels, 1)
        styles = _per_level(kw.get("ls"), levels, "-")
        view_lo = np.array([g.x[0], g.y[0]])
        view_hi = np.array([g.x[-1], g.y[-1]])
        for geom in snap.hier.levels:
            if geom.level not in levels or geom.level > g.level:
                continue
            segs = patch_outline_segments(geom, view_lo, view_hi)
            lc = LineCollection(segs, colors=colors[geom.level], linewidths=widths[geom.level],
                                linestyles=styles[geom.level], alpha=0.4)
            lc.set_label(f"L{geom.level} patches")
            ax.add_collection(lc)

    ax.set_aspect(kw.get("aspect", "equal"))
    ax.set_title(kw.get("title", default_title(field)))
    ax.set_xlabel(kw.get("xlabel", "x"))
    ax.set_ylabel(kw.get("ylabel", "y"))
    if kw.get("xlim") is not None:
        ax.set_xlim(kw["xlim"])
    if kw.get("ylim") is not None:
        ax.set_ylim(kw["ylim"])
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.08)
    fig.colorbar(im, ax=ax, cax=cax, extend=extend)
    return fig, ax


# DESIGN-NOTE: master draws one line per patch labelled L<l>P<i>; here one line per level
# (labelled L<l>, NaN between patches) - same picture, no per-patch Python loop - and the 1D
# plot returns (fig, ax) like the 2D one (master returned None).
def _plot1d(field, kw, qty):
    snap = field.snapshot()
    present = [g.level for g in snap.hier.levels]
    levels = _levels(kw.get("levels"), present, default_all=False)
    fig, ax = _figure(kw.get("ax"))
    for geom, ld in zip(snap.hier.levels, snap.levels):
        if geom.level not in levels or geom.npatch == 0:
            continue
        # one line per level; NaN between patches breaks the line at patch ends
        values = np.asarray(ld.buffer("value"), dtype=np.float64)
        x = geom.node_coords(0)
        cut = geom.offset[1:-1]
        ax.plot(np.insert(x, cut, np.nan), np.insert(values, cut, np.nan),
                label=f"L{geom.level}", marker=kw.get("marker", ""), ls=kw.get("ls", "--"),
                color=kw.get("color", "k"))
    ax.set_title(kw.get("title", default_title(field)))
    ax.set_xlabel(kw.get("xlabel", "x"))
    ax.set_ylabel(kw.get("ylabel", qty if qty is not None else (field.quantity or "")))
    if kw.get("xlim") is not None:
        ax.set_xlim(kw["xlim"])
    if kw.get("ylim") is not None:
        ax.set_ylim(kw["ylim"])
    return fig, ax


def plot(field, **kwargs):
    """Draw a 1D or 2D field with master's keywords; returns ``(fig, ax)``.

    Master's keywords: qty (component: "x", "Bx", "xx"...; required for vector/tensor
    fields), time (needed for a multi-time field), levels (1D default (0,); 2D default all:
    ``levels=[1]`` shows only the nodes level 1 holds, the rest is masked), ax, cmap, vmin,
    vmax, plot_patches (2D outlines, one LineCollection per level, only patches inside the
    view), patchcolors, lw, ls, aspect, title, xlabel, ylabel, xlim, ylim (2D: the zoom;
    only the patches under it are read), legend, filename (saved with ``dpi``, default 200),
    dpi, marker, color (1D). A component field accepts its parent's name as qty
    (``B.x.plot(qty="x")``). Default title: empty, or the quantity and its stencil for a
    derived field (``"J_z (stencil: node)"``).
    New keywords: robust (percentile colour limits), symmetric, norm (a matplotlib Normalize,
    instead of vmin/vmax), max_level (finest level drawn), fill ("cubic" or "nearest", how
    finest() fills nodes no level holds).
    Reads only what is drawn. Raises TypeError for an unknown keyword, ValueError for a 3D
    field (take a cut first), a vector field without qty, or norm with vmin/vmax.
    Master returned None in 1D; this returns ``(fig, ax)`` in 1D and 2D.
    """
    _check_keywords(kwargs)
    if field.ndim == 3:
        raise ValueError(
            "plot() draws 1D and 2D fields; for a 3D field take a cut first, e.g. "
            "g = field.finest(box=dict(z=(z0, z0))) and imshow a slice of g.data, or a line "
            "probe field(y=y0, z=z0)"
        )
    field = _one_time(field, kwargs.get("time"))
    comp = _component(field, kwargs.get("qty"))
    if field.ndim == 1:
        fig, ax = _plot1d(comp, kwargs, kwargs.get("qty"))
    else:
        fig, ax = _plot2d(comp, kwargs)
    if kwargs.get("legend") is not None:
        ax.legend()
    if kwargs.get("filename") is not None:
        fig.savefig(kwargs["filename"], dpi=kwargs.get("dpi", 200))
    return fig, ax


def install():
    """Make ``plot`` a method of every AMR field (``AMRField.plot = plot``)."""
    from ..field import AMRField

    AMRField.plot = plot


__all__ = ["MASTER_KEYWORDS", "NEW_KEYWORDS", "default_title", "install",
           "patch_outline_segments", "plot"]
