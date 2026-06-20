#
#
#


import numpy as np


def plot_field_data(fd, **kwargs):
    """Plot a single FieldData object."""
    from . import get_fig_ax

    fig, ax = get_fig_ax(**kwargs)

    if fd.ndim == 1:
        x, data = _strip_ghosts_1d(fd, fd.box)
        plot_kwargs = {
            k: kwargs[k] for k in ("label", "color", "ls", "marker") if k in kwargs
        }
        ax.plot(x, data, **plot_kwargs)
        ax.set_xlabel(kwargs.get("xlabel", "x"))
        ax.set_ylabel(kwargs.get("ylabel", fd.field_name))

    elif fd.ndim == 2:
        x, y, data = _pcolormesh_coords(fd, fd.box)
        im = ax.pcolormesh(
            x,
            y,
            data.T,
            cmap=kwargs.get("cmap", "Spectral_r"),
            vmin=kwargs.get("vmin", None),
            vmax=kwargs.get("vmax", None),
        )
        ax.set_aspect(kwargs.get("aspect", "equal"))
        ax.set_xlabel(kwargs.get("xlabel", "x"))
        ax.set_ylabel(kwargs.get("ylabel", "y"))
        _add_colorbar(ax, im)

    kwargs["ax"] = ax
    _finalize_ax(fig, **kwargs)
    return fig, ax


def plot_patches(hier, save=False):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(10, 3))
    for ilvl, lvl in hier.levels(0.0).items():
        lvl_offset = ilvl * 0.1
        for patch in lvl.patches:
            dx = patch.dl[0]
            x0 = patch.box.lower * dx
            x1 = patch.box.upper * dx
            xcells = np.arange(x0, x1 + dx, dx)
            y = lvl_offset + np.zeros_like(xcells)
            ax.plot(xcells, y, marker=".")
    if save:
        fig.savefig("hierarchy.png")


def plot_2d_patches(hier, ilvl, collections, **kwargs):
    import matplotlib.pyplot as plt
    from matplotlib.collections import PatchCollection

    from ....core import box as boxm
    from . import box_to_Rectangle

    if isinstance(collections, list) and all(
        isinstance(el, boxm.Box) for el in collections
    ):
        collections = [{"boxes": collections}]

    level_domain_box = hier.level_domain_box(ilvl)
    mi, ma = level_domain_box.lower.min(), level_domain_box.upper.max()

    fig, ax = kwargs.get("subplot", plt.subplots(figsize=(6, 6)))

    for collection in collections:
        facecolor = collection.get("facecolor", "none")
        edgecolor = collection.get("edgecolor", "purple")
        alpha = collection.get("alpha", 1)
        rects = [box_to_Rectangle(box) for box in collection["boxes"]]
        ax.add_collection(
            PatchCollection(
                rects, facecolor=facecolor, alpha=alpha, edgecolor=edgecolor
            )
        )

    if "title" in kwargs:
        from textwrap import wrap

        xfigsize = int(fig.get_size_inches()[0] * 10)
        ax.set_title("\n".join(wrap(kwargs["title"], xfigsize)))

    major_ticks = np.arange(mi - 5, ma + 5 + 5, 5)
    ax.set_xticks(major_ticks)
    ax.set_yticks(major_ticks)
    minor_ticks = np.arange(mi - 5, ma + 5 + 5, 1)
    ax.set_xticks(minor_ticks, minor=True)
    ax.set_yticks(minor_ticks, minor=True)

    ax.grid(which="both")
    return fig


def plot1d(hier, **kwargs):
    from . import get_fig_ax

    usr_lvls = kwargs.get("levels", (0,))
    qty = kwargs.get("qty", None)
    time = kwargs.get("time", hier.times()[0])

    fig, ax = get_fig_ax(**kwargs)

    if qty is None:
        if len(hier.quantities()) != 1:
            raise ValueError(
                "multiple quantities in patch, please specify a quantity in "
                + " ".join(hier.quantities())
            )
        qty = hier.quantities()[0]

    for lvl_nbr, level in hier.levels(time).items():
        if lvl_nbr not in usr_lvls:
            continue
        for ip, patch in enumerate(level.patches):
            x, val = _strip_ghosts_1d(patch[qty], patch.box)
            ax.plot(
                x,
                val,
                label=f"L{lvl_nbr}P{ip}",
                marker=kwargs.get("marker", ""),
                ls=kwargs.get("ls", "--"),
                color=kwargs.get("color", "k"),
            )

    ax.set_xlabel(kwargs.get("xlabel", "x"))
    ax.set_ylabel(kwargs.get("ylabel", qty))

    kwargs["ax"] = ax
    _finalize_ax(fig, **kwargs)
    return fig, ax


def plot2d(hier, **kwargs):
    from matplotlib.patches import Rectangle

    from . import get_fig_ax

    time = kwargs.get("time", hier._default_time())
    usr_lvls = kwargs.get("levels", hier.levelNbrs(time))
    default_qty = hier.quantities()[0] if len(hier.quantities()) == 1 else None
    qty = kwargs.get("qty", default_qty)

    fig, ax = get_fig_ax(**kwargs)

    glob_min = hier.global_min(qty)
    glob_max = hier.global_max(qty)

    patchcolors = kwargs.get("patchcolors", {ilvl: "k" for ilvl in usr_lvls})
    if not isinstance(patchcolors, dict):
        patchcolors = dict(zip(usr_lvls, patchcolors))

    linewidths = kwargs.get("lw", {ilvl: 1 for ilvl in usr_lvls})
    if not isinstance(linewidths, dict):
        linewidths = dict(zip(usr_lvls, linewidths))

    linestyles = kwargs.get("ls", {ilvl: "-" for ilvl in usr_lvls})
    if not isinstance(linestyles, dict):
        linestyles = dict(zip(usr_lvls, linestyles))

    for lvl_nbr, _ in hier.levels(time).items():
        if lvl_nbr not in usr_lvls:
            continue
        for patch in hier.level(lvl_nbr, time).patches:
            pdat = patch[qty]
            x, y, data = _pcolormesh_coords(pdat, patch.box)
            dx, dy = pdat.dl
            im = ax.pcolormesh(
                x,
                y,
                data.T,
                cmap=kwargs.get("cmap", "Spectral_r"),
                vmin=kwargs.get("vmin", glob_min - 1e-6),
                vmax=kwargs.get("vmax", glob_max + 1e-6),
            )
            if kwargs.get("plot_patches", False):
                r = Rectangle(
                    (patch.box.lower[0] * dx, patch.box.lower[1] * dy),
                    patch.box.shape[0] * dx,
                    patch.box.shape[1] * dy,
                    fc="none",
                    ec=patchcolors[lvl_nbr],
                    alpha=0.4,
                    lw=linewidths[lvl_nbr],
                    ls=linestyles[lvl_nbr],
                )
                ax.add_patch(r)

    ax.set_aspect(kwargs.get("aspect", "equal"))
    ax.set_xlabel(kwargs.get("xlabel", "x"))
    ax.set_ylabel(kwargs.get("ylabel", "y"))
    _add_colorbar(ax, im)

    kwargs["ax"] = ax
    _finalize_ax(fig, **kwargs)
    return fig, ax


def plot(hier, **kwargs):
    if hier.ndim == 1:
        return plot1d(hier, **kwargs)
    elif hier.ndim == 2:
        return plot2d(hier, **kwargs)
    raise ValueError("3d not supported for plotting")


def _strip_ghosts_1d(pdat, box):
    ng = pdat.ghosts_nbr[0]
    if ng > 0:
        return pdat.x[ng:-ng], pdat[box]
    return pdat.x, pdat.dataset[:]


def _pcolormesh_coords(pdat, box):
    """Strip ghosts and return (x, y, data) with edge coordinates for pcolormesh."""
    ng = pdat.ghosts_nbr
    data = pdat[box] if np.any(ng != 0) else pdat.dataset[:]
    sx = slice(ng[0], -ng[0] if ng[0] else None)
    sy = slice(ng[1], -ng[1] if ng[1] else None)
    x = np.copy(pdat.x[sx])
    y = np.copy(pdat.y[sy])
    return x, y, data


def _add_colorbar(ax, im):
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.08)
    plt.colorbar(im, ax=ax, cax=cax)


def _finalize_ax(fig, **kwargs):
    ax = kwargs["ax"]
    ax.set_title(kwargs.get("title", ""))
    if "xlim" in kwargs:
        ax.set_xlim(kwargs["xlim"])
    if "ylim" in kwargs:
        ax.set_ylim(kwargs["ylim"])
    if kwargs.get("legend") is not None:
        ax.legend()
    if "filename" in kwargs:
        fig.savefig(kwargs["filename"], dpi=kwargs.get("dpi", 200))
