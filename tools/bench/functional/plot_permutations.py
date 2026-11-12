#
# plots the permutor results, after:
#   python3 tools/python3/bench.py print_summary -i .phare_bench/<datetime>.tar.gz > perms
#   python3 -c "from tools.bench.functional.permutor import sort_summaries as sort; sort('perms')" # optional
# then:
#   python3 tools/bench/functional/plot_permutations.py -i perms
#
# one set of figures per combination of the permuted keys other than -x and -l
# perms_time.png: per level values over the run, from the permutor's summary samples
#

import argparse
import sys
from pathlib import Path

import numpy as np

from tools.python3.bench import parse_stats_line


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", default="perms", help="print_summary output")
    parser.add_argument(
        "-x", "--x-axis", default="tag_buffer", help="permuted key on the x axis"
    )
    parser.add_argument(
        "-l", "--line", default="tagging_threshold", help="permuted key per line"
    )
    parser.add_argument("-o", "--outdir", default=".", help="where to save the figures")
    args = parser.parse_args()

    with open(args.input) as f:
        records = [r for r in (parse_stats_line(line) for line in f) if r]
    if not records:
        sys.exit(f"no print_summary lines found in {args.input}")
    for flag, key in (("--x-axis", args.x_axis), ("--line", args.line)):
        # list values like cells split the groups instead
        if any(isinstance(r["kwargs"].get(key), list) for r in records):
            sys.exit(f"{flag} needs a scalar permuted key, {key} has list values")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    for name, group in groups(records, args.x_axis, args.line).items():
        plot_group(group, args.x_axis, args.line, outdir, name)
    print_fastest(records)


def groups(records, x, line):
    """records by their values of the permuted keys other than x and line"""
    keys = [
        k
        for k in records[0]["kwargs"]
        if k not in (x, line) and len({repr(r["kwargs"].get(k)) for r in records}) > 1
    ]
    grouped = {}
    for r in records:
        name = "_".join(f"{k}={_compact(r['kwargs'].get(k))}" for k in keys)
        grouped.setdefault(name, []).append(r)
    return grouped


def plot_group(records, x, line, outdir, name):
    import matplotlib.pyplot as plt

    suffix = f"_{name}" if name else ""
    xs = sorted({r["kwargs"][x] for r in records})
    ls = sorted({r["kwargs"][line] for r in records})

    # run wide values against x, one line per value of `line`
    metrics = [
        ("ppp_ns", "push time per particle (ns)"),
        ("time_s", "average L0 advance (s)"),
        ("mem_mb", "max memory (MB)"),
        ("particles", "particles"),
    ]
    fig, axes = plt.subplots(1, len(metrics), figsize=(4 * len(metrics), 3.5))
    fig.suptitle(name or "permutations")
    for ax, (key, label) in zip(axes, metrics):
        for lv in ls:
            pts = _points(records, x, line, lv, lambda r, key=key: r[key])
            ax.plot(*pts, marker="o", label=f"{line}={lv}")
        _decorate(ax, x, label, xs)
    axes[-1].legend(fontsize=7, bbox_to_anchor=(1.02, 1), loc="upper left")
    _save(fig, outdir / f"perms_runs{suffix}.png")

    # push time heatmap, x by line
    mat = np.full((len(ls), len(xs)), np.nan)
    for r in records:
        if r["ppp_ns"] is not None:
            mat[ls.index(r["kwargs"][line]), xs.index(r["kwargs"][x])] = r["ppp_ns"]
    fig, ax = plt.subplots(figsize=(1.2 * len(xs) + 2, 0.6 * len(ls) + 2))
    im = ax.imshow(mat, aspect="auto", cmap="RdYlGn_r")
    ax.set_xticks(range(len(xs)), [str(v) for v in xs])
    ax.set_yticks(range(len(ls)), [str(v) for v in ls])
    ax.set(xlabel=x, ylabel=line, title=f"push time per particle (ns) {name}")
    for (i, j), v in np.ndenumerate(mat):
        if not np.isnan(v):
            ax.text(j, i, f"{v:.1f}", ha="center", va="center", fontsize=7)
    fig.colorbar(im, ax=ax)
    _save(fig, outdir / f"perms_ppp_heatmap{suffix}.png")

    # per level values, one column per level
    layout = [  # (label, [(key, linestyle), ...])
        (
            "% of level refined",
            [("refined_pc", "-")],
        ),  # level i: % covered by level i + 1
        ("ghost / domain cells", [("cell_ratio", "-")]),
        ("cell imbalance (max/avg)", [("cell_imbalance", "-")]),
        ("particle imbalance (max/avg)", [("particle_imbalance", "-")]),
    ]
    _level_figure(
        records,
        x,
        line,
        xs,
        ls,
        layout,
        f"per level {name}",
        outdir / f"perms_levels{suffix}.png",
    )

    # patch ghost exchange: cells filled from same level neighbours on other ranks
    mpi = [
        ("remote patch ghost cells", [("remote_ghosts", "-")]),
        ("% patch ghosts remote", [("remote_pc", "-")]),
        (
            "remote ranks (- max, -- avg)",
            [("remote_ranks_max", "-"), ("remote_ranks_avg", "--")],
        ),
    ]
    _level_figure(
        records,
        x,
        line,
        xs,
        ls,
        mpi,
        f"MPI patch ghost exchange per level {name}",
        outdir / f"perms_mpi{suffix}.png",
    )

    # the same per level values over the run, from the permutor's summary samples
    over_time = [
        ("cells", "cells"),
        ("cell_imbalance", "cell imbalance (max/avg)"),
        ("particle_imbalance", "particle imbalance (max/avg)"),
        ("remote_ghosts", "remote patch ghost cells"),
        ("remote_ranks_max", "max remote ranks"),
    ]
    _time_figure(
        records,
        x,
        line,
        xs,
        ls,
        over_time,
        f"per level over time {name}",
        outdir / f"perms_time{suffix}.png",
    )


def print_fastest(records, n=10):
    timed = sorted((r for r in records if r["ppp_ns"]), key=lambda r: r["ppp_ns"])
    print(f"\nTop {min(n, len(timed))} fastest (push time per particle):")
    for r in timed[:n]:
        print(
            f"  {r['ppp_ns']:.1f} ns  MEM={r['mem_mb']}MB  N={r['particles']:,}  {r['kwargs']}"
        )


def _level_figure(records, x, line, xs, ls, rows, title, path):
    """rows of (label, [(key, linestyle)]), one column per level, one colour per line value"""
    import matplotlib.pyplot as plt

    n_levels = max(len(r["levels"]) for r in records)
    fig, axes = plt.subplots(
        len(rows), n_levels, figsize=(4 * n_levels, 2.5 * len(rows)), squeeze=False
    )
    fig.suptitle(title)
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    for ilvl in range(n_levels):
        for (label, keys), ax in zip(rows, axes[:, ilvl]):
            for i, lv in enumerate(ls):
                for key, linestyle in keys:
                    pts = _points(
                        records,
                        x,
                        line,
                        lv,
                        lambda r, key=key, ilvl=ilvl: _level(r["levels"], ilvl, key),
                    )
                    ax.plot(
                        *pts,
                        marker="o",
                        linestyle=linestyle,
                        color=colors[i % len(colors)],
                        label=f"{line}={lv}" if linestyle == "-" else None,
                    )
            _decorate(ax, x, label, xs)
            ax.set_visible(any(len(plotted.get_xdata()) for plotted in ax.lines))
        # on the first panel with data, e.g. not % refined on the finest level
        next((ax for ax in axes[:, ilvl] if ax.get_visible()), axes[0, ilvl]).set_title(
            f"level {ilvl}"
        )
    # on the first panel of the last level with data, e.g. not % refined on the finest
    legend_ax = next((ax for ax in axes[:, -1] if ax.get_visible()), axes[0, -1])
    legend_ax.legend(fontsize=7, bbox_to_anchor=(1.02, 1), loc="upper left")
    _save(fig, path)


def _time_figure(records, x, line, xs, ls, rows, title, path):
    """rows of (key, label) against time, one column per level, one line per run:
    colour per value of `line`, line style per value of `x`"""
    import matplotlib.pyplot as plt

    n_levels = max(len(s["levels"]) for r in records for s in r["samples"])
    fig, axes = plt.subplots(
        len(rows), n_levels, figsize=(4 * n_levels, 2.5 * len(rows)), squeeze=False
    )
    fig.suptitle(title)
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    styles = ["-", "--", ":", "-."]
    for ilvl in range(n_levels):
        for (key, label), ax in zip(rows, axes[:, ilvl]):
            for r in sorted(records, key=lambda r: (r["kwargs"][x], r["kwargs"][line])):
                pts = [
                    (s["time"], _level(s["levels"], ilvl, key)) for s in r["samples"]
                ]
                pts = [(t, v) for t, v in pts if v is not None]
                if not pts:
                    continue
                xi, li = xs.index(r["kwargs"][x]), ls.index(r["kwargs"][line])
                ax.plot(
                    *zip(*pts),
                    marker=".",
                    color=colors[li % len(colors)],
                    linestyle=styles[xi % len(styles)],
                    label=f"{x}={xs[xi]} {line}={ls[li]}",
                )
            ax.set(xlabel="time", ylabel=label)
            ax.grid(True, alpha=0.3)
            ax.set_visible(bool(ax.lines))
        next((ax for ax in axes[:, ilvl] if ax.get_visible()), axes[0, ilvl]).set_title(
            f"level {ilvl}"
        )
    legend_ax = next((ax for ax in axes[:, -1] if ax.get_visible()), axes[0, -1])
    legend_ax.legend(fontsize=6, bbox_to_anchor=(1.02, 1), loc="upper left")
    _save(fig, path)


def _points(records, x, line, line_value, value_of):
    pts = sorted(
        (r["kwargs"][x], value_of(r))
        for r in records
        if r["kwargs"][line] == line_value
    )
    pts = [(xv, yv) for xv, yv in pts if yv is not None]
    return [p[0] for p in pts], [p[1] for p in pts]


def _level(levels, ilvl, key):  # None when the level or the value is missing
    return levels[ilvl].get(key) if ilvl < len(levels) else None


def _decorate(ax, x, label, xs):
    ax.set(xlabel=x, ylabel=label)
    ax.set_xticks(xs)
    ax.grid(True, alpha=0.3)


def _save(fig, path):
    import matplotlib.pyplot as plt

    fig.tight_layout()
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    print("saved", path)


def _compact(value):  # for file names
    return (
        str(value).replace(" ", "").replace("[", "").replace("]", "").replace(",", "x")
    )


if __name__ == "__main__":
    main()
