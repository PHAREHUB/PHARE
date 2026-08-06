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
