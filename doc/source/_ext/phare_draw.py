"""Shared look for the drawings of the developer guide (matplotlib, drawn at build time).

Use from a ``.. plot::`` block::

    .. plot::
       :include-source: False
       :caption: one line: what to look at

       from phare_draw import *
       fig, (a1, a2) = strip(2, titles=["Coarse cell", "Fine cells"])
       grid(a1, (0, 0, 4, 4)); yee(a1, 1, 1)
       finish(fig, [a1, a2])      # always last: sizes every panel to its content (no empty margins)

Conventions (all lengths are in *cells*; a cell i spans [i, i+1)):
  primal = integer position (a node, on a cell boundary), dual = half-integer position (cell centre).
  Yee layout in 2D, as in src/core/data/grid/impl/yee/gridlayout_hybrid_yee.hpp:
    Bx (primal x, dual y)  By (dual x, primal y)  Bz (dual, dual)
    Ex (dual, primal)      Ey (primal, dual)      Ez (primal, primal)     density/moments: (primal, primal)
Rules: width <= 6.5 in, panels stacked vertically, fontsize >= 11. Figures have an opaque white
background, so they read the same on a light or a dark page.
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyArrowPatch, Rectangle, FancyBboxPatch  # noqa: E402

__all__ = [
    "C", "strip", "grid", "patch", "level", "yee", "mark", "arrow", "particles", "text",
    "brace", "boxes", "frame", "finish", "WIDTH", "legend", "ruler", "yee_pair",
]

WIDTH = 6.5  # inches, maximum
HALO = dict(fc="white", ec="none", pad=0.8, alpha=0.9)  # label background: text never sits on a line

# Palette: mid-saturation colours that keep contrast on white.
C = dict(
    ink="#212529", grey="#868e96", light="#dee2e6", faint="#f1f3f5",
    B="#d6336c", E="#1971c2", J="#e67700", ion="#2f9e44", rho="#7048e8",
    ghost="#ffe8cc", ghost_edge="#e8590c", patch="#e7f5ff", patch_edge="#1971c2",
    L0="#868e96", L1="#1971c2", L2="#d9480f", ok="#2f9e44", bad="#e03131", hi="#fff3bf",
)

plt.rcParams.update({
    "font.size": 11, "axes.titlesize": 12, "axes.titleweight": "bold", "axes.titlelocation": "left",
    "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white",
    "axes.edgecolor": C["grey"], "text.color": C["ink"], "axes.labelcolor": C["ink"],
    "font.family": "DejaVu Sans", "mathtext.default": "regular", "lines.solid_capstyle": "round",
})


def strip(n=1, titles=None, height=2.6, width=WIDTH, aspect=True, hspace=0.35, phone=False):
    """n panels stacked vertically, titled "1 · ...", "2 · ..." (a step-by-step strip).

    phone=True makes the figure 3.4 in wide (340 px at 100 dpi): it is seen whole on a phone. Use it when a panel
    has few elements; dense panels keep the default 6.5 in (they scroll sideways on a phone).

    Returns (fig, axes): axes is a list even when n == 1. Axes are 'off' and (by default) equal aspect.
    """
    if phone:
        width = 3.4
    fig, axes = plt.subplots(n, 1, figsize=(min(width, WIDTH), height * n), squeeze=False)
    axes = list(axes[:, 0])
    for i, ax in enumerate(axes):
        ax.set_axis_off()
        if aspect:
            ax.set_aspect("equal")
        if titles:
            t = titles[i]
            ax._phare_title = f"{i + 1} · {t}" if n > 1 else t
    fig.subplots_adjust(left=0.02, right=0.98, top=0.93 if n == 1 else 0.96, bottom=0.02, hspace=hspace)
    return fig, axes


def grid(ax, box, color=None, lw=0.8, fill=None, zorder=1):
    """Cell grid over box = (x0, y0, x1, y1) in cell units."""
    x0, y0, x1, y1 = box
    color = color or C["light"]
    if fill:
        ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fc=fill, ec="none", zorder=zorder - 1))
    for x in range(int(x0), int(x1) + 1):
        ax.plot([x, x], [y0, y1], color=color, lw=lw, zorder=zorder)
    for y in range(int(y0), int(y1) + 1):
        ax.plot([x0, x1], [y, y], color=color, lw=lw, zorder=zorder)
    ax.update_datalim([(x0, y0), (x1, y1)])
    ax.autoscale_view()


def patch(ax, box, ghosts=0, label=None, color=None, ghost_color=None, lw=2.2, dl=1.0):
    """A patch: box = (x0, y0, x1, y1) in cell edges, interior filled light blue, `ghosts` cells of ghost layer around it
    (light orange tint, dashed outline). Returns the outer box including ghosts."""
    x0, y0, x1, y1 = box
    ec = color or C["patch_edge"]
    if ghosts:
        g = ghosts
        ax.add_patch(Rectangle((x0 - g, y0 - g), x1 - x0 + 2 * g, y1 - y0 + 2 * g, fc=ghost_color or C["ghost_edge"],
                               ec=C["ghost_edge"], lw=1.4, ls="--", zorder=0.5, alpha=0.15))
        ax.add_patch(Rectangle((x0 - g, y0 - g), x1 - x0 + 2 * g, y1 - y0 + 2 * g, fc="none",
                               ec=C["ghost_edge"], lw=1.4, ls="--", zorder=0.6))
    ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fc=C["patch"], ec="none", zorder=0.8))
    ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fc="none", ec=ec, lw=lw, zorder=1.5))  # above grid lines
    if label:
        ax.text((x0 + x1) / 2, y1 + (ghosts + 0.25), label, ha="center", va="bottom", color=ec, fontweight="bold")
    g = ghosts
    ax.update_datalim([(x0 - g, y0 - g), (x1 + g, y1 + g)])
    ax.autoscale_view()
    return (x0 - g, y0 - g, x1 + g, y1 + g)


def level(ax, box, ratio_label=None, color=None, lw=2.5, ls="-", ratio=1):
    """Outline of a refined region (a level or a group of patches), box in cell edges of *this panel's* grid.
    ratio=2 (or 4) also draws the fine cells inside the box, so the refinement is visible."""
    x0, y0, x1, y1 = box
    if ratio > 1:
        n = ratio
        step = 1.0 / n
        k = 0
        while x0 + k * step <= x1 + 1e-9:
            ax.plot([x0 + k * step] * 2, [y0, y1], color=color or C["L1"], lw=0.6, alpha=0.55, zorder=2.5)
            k += 1
        k = 0
        while y0 + k * step <= y1 + 1e-9:
            ax.plot([x0, x1], [y0 + k * step] * 2, color=color or C["L1"], lw=0.6, alpha=0.55, zorder=2.5)
            k += 1
    ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fc="none", ec=color or C["L1"], lw=lw, ls=ls, zorder=3))
    if ratio_label:
        ax.text(x0, y1 + 0.2, ratio_label, color=color or C["L1"], fontweight="bold", va="bottom", bbox=HALO)


def mark(ax, kind, x, y, color=None, size=9, label=None, zorder=5):
    """Value markers. kind: 'primal' (filled circle), 'dual' (open square), 'cell' (cross), 'node' (small dot)."""
    color = color or C["ink"]
    if kind == "primal":
        ax.plot(x, y, "o", ms=size, color=color, mec="white", mew=0.8, zorder=zorder)
    elif kind == "dual":
        ax.plot(x, y, "s", ms=size, mfc="white", mec=color, mew=1.8, zorder=zorder)
    elif kind == "cell":
        ax.plot(x, y, "x", ms=size, color=color, mew=2, zorder=zorder)
    else:
        ax.plot(x, y, ".", ms=size / 2, color=color, zorder=zorder)
    if label:
        ax.text(x + 0.12, y + 0.12, label, color=color, fontsize=11, zorder=zorder)


def arrow(ax, p, q, color=None, lw=2, text=None, style="-|>", ms=14, ls="-", zorder=6, text_off=(0, 0.12),
          text_ha="center", **kw):
    """Arrow from p to q (data coordinates), optional label at the middle."""
    a = FancyArrowPatch(p, q, arrowstyle=style, mutation_scale=ms, lw=lw, color=color or C["ink"], ls=ls, zorder=zorder, **kw)
    ax.add_patch(a)
    if text:
        ax.text((p[0] + q[0]) / 2 + text_off[0], (p[1] + q[1]) / 2 + text_off[1], text, ha=text_ha,
                va="bottom" if text_off[1] else "center", color=color or C["ink"], bbox=HALO, zorder=zorder + 1)
    return a


# where each Yee component sits in 2D, relative to the cell corner (i, j), and which way it points
_YEE = {
    "Bx": (0.0, 0.5, (1, 0), "B"), "By": (0.5, 0.0, (0, 1), "B"), "Bz": (0.5, 0.5, None, "B"),
    "Ex": (0.5, 0.0, (1, 0), "E"), "Ey": (0.0, 0.5, (0, 1), "E"), "Ez": (0.0, 0.0, None, "E"),
}


# label position inside the cell (relative to its corner), chosen so no two labels of one cell touch
_LABEL = {
    "Bx": (0.05, 0.57, "left", "bottom"), "By": (0.58, 0.10, "left", "bottom"), "Bz": (0.60, 0.57, "left", "bottom"),
    "Ex": (0.50, 0.08, "center", "bottom"), "Ey": (0.10, 0.50, "left", "center"), "Ez": (0.10, 0.08, "left", "bottom"),
}


def yee(ax, i, j, comps=("Bx", "By", "Bz", "Ex", "Ey", "Ez"), size=9, labels=True, fs=11):
    """Draw the Yee components of cell (i, j) at their PHARE places (see module docstring).
    B and E are arrows along their direction (Bx, Ex along x; By, Ey along y), Bz and Ez are dots (out of the plane).
    Arrows stay inside the cell; labels sit on the cell-interior side of their marker, on a white halo.
    Bx and Ey (and By and Ex) share a place: draw B and E in two cells with yee_pair()."""
    for name in comps:
        dx, dy, d, kind = _YEE[name]
        x, y = i + dx, j + dy
        col = C[kind]
        if d is None:
            ax.plot(x, y, "o", ms=size, mfc=col, mec="white", zorder=5)
            ax.plot(x, y, ".", ms=3, color="white", zorder=6)
        else:
            h = 0.17
            hx, hy = d[0] * h, d[1] * h
            # keep the arrow inside the cell: shift it inward if it sits on the cell boundary
            sx = 0.17 if (dx == 0.0 and d[0] == 1) else 0.0
            sy = 0.17 if (dy == 0.0 and d[1] == 1) else 0.0
            arrow(ax, (x + sx - hx, y + sy - hy), (x + sx + hx, y + sy + hy), color=col, lw=2.2, ms=12)
            # a little tick marks the exact position
            ax.plot(x, y, "o", ms=4, color=col, zorder=5)
        if labels:
            lx, ly, ha, va = _LABEL[name]
            ax.text(i + lx, j + ly, name, color=col, fontsize=fs, fontweight="bold", zorder=7, ha=ha, va=va, bbox=HALO)


def yee_pair(ax, origin=(0, 0), gap=1, comps_b=("Bx", "By", "Bz"), comps_e=("Ex", "Ey", "Ez")):
    """B and E on two copies of ONE cell, side by side (they share places, so one cell would overlap).
    Left: B on the faces and in the centre. Right: E on the edges and at the node. Returns the bounding box."""
    x, y = origin
    for k, (comps, name, kind) in enumerate(((comps_b, "B", "B"), (comps_e, "E", "E"))):
        ox = x + k * (1 + gap)
        grid(ax, (ox, y, ox + 1, y + 1), color=C["grey"], lw=1.2, fill=C["faint"])
        yee(ax, ox, y, comps=comps)
        ax.text(ox + 0.5, y + 1.12, name, ha="center", va="bottom", fontweight="bold", color=C[kind])
    return (x, y, x + 2 + gap, y + 1.4)


def legend(ax, entries, loc="center left", bbox=(1.0, 0.5), ncol=1, fs=11):
    """Legend outside the drawing, beside it. entries: list of (label, color, marker) with marker in
    'primal' (filled dot), 'dual' (open square), 'particle' (small dot), 'line', 'patch'."""
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    hs = []
    for label, color, m in entries:
        if m == "primal":
            hs.append(Line2D([], [], marker="o", ls="", color=color, ms=9, label=label))
        elif m == "dual":
            hs.append(Line2D([], [], marker="s", ls="", mfc="white", mec=color, mew=1.8, ms=9, label=label))
        elif m == "particle":
            hs.append(Line2D([], [], marker="o", ls="", color=color, ms=5, label=label))
        elif m == "line":
            hs.append(Line2D([], [], color=color, lw=2, label=label))
        else:
            hs.append(Patch(fc=color, ec=color, alpha=0.4, label=label))
    return ax.legend(handles=hs, loc=loc, bbox_to_anchor=bbox, frameon=False, fontsize=fs, ncol=ncol,
                     handletextpad=0.4, borderaxespad=0.2)


def ruler(ax, x0, x1, y=0.0, primal=True, dual=True, ticks=True, labels=True):
    """A 1D grid line from x0 to x1 (cell edges): primal values (filled dots) ON the nodes (cell boundaries), dual values (open
    squares) at the cell centres, ON the same line. labels: write the cell number under each centre."""
    ax.plot([x0, x1], [y, y], color=C["grey"], lw=1.4, zorder=1)
    ax.update_datalim([(x0, y - 0.8), (x1, y + 1.1)])  # room for a legend above and labels below
    for x in range(int(x0), int(x1) + 1):
        ax.plot([x, x], [y - 0.12, y + 0.12], color=C["grey"], lw=1.4, zorder=1)
        if primal:
            mark(ax, "primal", x, y, color=C["rho"])
    for x in range(int(x0), int(x1)):
        if dual:
            mark(ax, "dual", x + 0.5, y, color=C["B"])
        if labels:
            ax.text(x + 0.5, y - 0.3, str(x), ha="center", va="top", color=C["grey"], fontsize=11)


def particles(ax, xy, color=None, size=26, marker="o", zorder=5, **kw):
    """Particles as dots; xy is a list of (x, y) in cell units (use y=const*0.5 for a 1D picture)."""
    xs, ys = zip(*xy)
    ax.scatter(xs, ys, s=size, c=color or C["ion"], marker=marker, zorder=zorder, edgecolors="white", linewidths=0.6, **kw)


def text(ax, x, y, s, color=None, ha="center", va="center", fs=11, bold=False, halo=True, **kw):
    ax.text(x, y, s, color=color or C["ink"], ha=ha, va=va, fontsize=fs, fontweight="bold" if bold else "normal",
            bbox=HALO if halo else None, zorder=8, **kw)


def brace(ax, x0, x1, y, label, color=None, up=True, h=0.15):
    """Dimension marker: a bracket under (or over) x0..x1 at height y, with a label."""
    color = color or C["grey"]
    s = 1 if up else -1
    ax.plot([x0, x0, x1, x1], [y, y + s * h, y + s * h, y], color=color, lw=1.4)
    ax.text((x0 + x1) / 2, y + s * (h + 0.06), label, ha="center", va="bottom" if up else "top", color=color)


def frame(ax, pad=0.3):
    """Fit the axes around everything drawn, with a margin (call last)."""
    ax.autoscale_view()
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    ax.set_xlim(x0 - pad, x1 + pad)
    ax.set_ylim(y0 - pad, y1 + pad)


def finish(fig, axes, max_h=3.4, min_h=0.9, pad=0.25, title_h=0.38, gap=0.18):
    """Always call last. Fits each panel's view to what is drawn, then sets each panel's height from its
    data aspect (equal aspect) so there is no empty margin: panel height = width * (data height / data width),
    clamped to [min_h, max_h] inches. Panels with a title get room for it."""
    W = fig.get_size_inches()[0]
    inner_w = W * 0.96
    hs = []
    for ax in axes:
        frame(ax, pad)
        x0, x1 = ax.get_xlim()
        y0, y1 = ax.get_ylim()
        h = inner_w * (y1 - y0) / max(x1 - x0, 1e-9)
        hs.append(min(max(h, min_h), max_h))
    titled = [bool(getattr(ax, "_phare_title", "")) for ax in axes]
    total = sum(hs) + sum(title_h if t else 0 for t in titled) + gap * (len(axes) - 1) + 0.1
    fig.set_size_inches(W, total)
    y = total - 0.05
    for ax, h, t in zip(axes, hs, titled):
        if t:
            fig.text(0.02, (y - 0.03) / total, ax._phare_title, ha="left", va="top", fontsize=11, fontweight="bold")
        y -= (title_h if t else 0)
        ax.set_position([0.02, (y - h) / total, 0.96, h / total])
        y -= h + gap
    return fig


def boxes(ax, items, arrows=(), w=2.6, h=0.9, gap=0.7, color=None, cols=1, fs=11):
    """A column (or `cols` columns) of labelled rounded boxes, with arrows between box indices.
    items: list of str ("title" or "title\\nsmall text"). arrows: list of (i, j) or (i, j, "label").
    For navigation chains prefer the `.. flow::` directive (clickable)."""
    pos = []
    n = len(items)
    per = -(-n // cols)
    for k, s in enumerate(items):
        c, r = divmod(k, per)
        x, y = c * (w + 0.9), -r * (h + gap)
        pos.append((x, y))
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.15", fc=C["faint"],
                                    ec=color or C["grey"], lw=1.6, zorder=2))
        lines = s.split("\n")
        ax.text(x + w / 2, y + h / 2 + (0.14 if len(lines) > 1 else 0), lines[0], ha="center", va="center",
                fontweight="bold", fontsize=fs, zorder=3)
        if len(lines) > 1:
            ax.text(x + w / 2, y + h / 2 - 0.2, "\n".join(lines[1:]), ha="center", va="center", fontsize=fs - 1,
                    color=C["grey"], zorder=3)
    for a in arrows:
        i, j = a[0], a[1]
        (xi, yi), (xj, yj) = pos[i], pos[j]
        p = (xi + w / 2, yi) if yj < yi else (xi + w, yi + h / 2)
        q = (xj + w / 2, yj + h) if yj < yi else (xj, yj + h / 2)
        arrow(ax, p, q, color=C["grey"], lw=1.6, text=a[2] if len(a) > 2 else None,
              text_off=(0.12, 0.0), text_ha="left")
    ax.autoscale_view()
    ax.set_aspect("equal")
    return pos
