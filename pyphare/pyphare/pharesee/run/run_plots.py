import os

import matplotlib.pyplot as plt
import numpy as np

from .run import Run


def finest_field_plot(run_path, qty, **kwargs):
    """
    plot the given quantity (qty) at 'run_path' with only the finest data

    * run_path : the path of the run
    * qty : ['Bx', 'By', 'Bz', 'Ex', 'Ey', 'Ez', 'Fx', 'Fy', 'Fz', 'Vx', 'Vy', 'Vz', 'rho']

    kwargs:
    * ax : the handle for the fig axes
    * time : time (that should be in the time_stamps of the run)
    * interp : the type of interpolation for the interpolator
    * draw_style : steps-mid ou default... only for 1d
    * title : (str) title of the plot
    * xlabel, ylabel
    * xlim, ylim
    * filename : (str) if exists, save plot to figure under that name

    return value : fig,ax
    """

    from mpl_toolkits.axes_grid1 import make_axes_locatable

    from pyphare.pharesee.hierarchy.fromh5 import get_times_from_h5

    r = Run(run_path)

    time = kwargs.get("time", None)
    dim = r.GetDl("finest", time).shape[0]
    interp = kwargs.get("interp", "nearest")

    if qty in ["Bx", "By", "Bz"]:
        file = os.path.join(run_path, "EM_B.h5")
        if time is None:
            times = get_times_from_h5(file)
            time = times[0]
        grid = r.GetB(time, all_primal=False).finest(time, qty, interp)
    elif qty in ["Ex", "Ey", "Ez"]:
        file = os.path.join(run_path, "EM_E.h5")
        if time is None:
            times = get_times_from_h5(file)
            time = times[0]
        grid = r.GetE(time, all_primal=False).finest(time, qty, interp)
    elif qty in ["Vx", "Vy", "Vz"]:
        file = os.path.join(run_path, "ions_bulkVelocity.h5")
        if time is None:
            times = get_times_from_h5(file)
            time = times[0]
        grid = r.GetVi(time).finest(time, qty[1:], interp)
    elif qty == "rho":
        file = os.path.join(run_path, "ions_charge_density.h5")
        if time is None:
            times = get_times_from_h5(file)
            time = times[0]
        grid = r.GetNi(time).finest(time, interp=interp)
    elif qty in ("Jx", "Jy", "Jz"):
        file = os.path.join(run_path, "EM_B.h5")
        if time is None:
            times = get_times_from_h5(file)
            time = times[0]
        grid = r.GetJ(time, all_primal=False).finest(time, qty, interp)
    else:
        # ___ TODO : should also include the files for a given population
        raise ValueError(
            "qty should be in ['Bx', 'By', 'Bz', 'Ex', 'Ey', 'Ez', 'Fx', 'Fy', 'Fz', 'Vx', 'Vy', 'Vz', 'rho']"
        )

    if "ax" not in kwargs:
        fig, ax = plt.subplots()
    else:
        ax = kwargs["ax"]
        fig = ax.figure

    if dim == 1:
        drawstyle = kwargs.get("drawstyle", "steps-mid")

        ax.plot(grid.x, grid[:], drawstyle=drawstyle)
    elif dim == 2:
        x = np.copy(grid.x)
        y = np.copy(grid.y)
        dx = x[1] - x[0]
        dy = y[1] - y[0]

        # pcolormesh considers DATA_ij to be the center of the pixel
        # and X,Y are the corners so XY need to be made 1 value larger
        # and shifted around DATA_ij
        x -= dx / 2
        x = np.append(x, x[-1] + dx)
        y -= dy / 2
        y = np.append(y, y[-1] + dy)

        DATA = grid[:].T  # pcolormesh takes (y, x)

        vmin = kwargs.get("vmin", np.nanmin(DATA))
        vmax = kwargs.get("vmax", np.nanmax(DATA))
        cmap = kwargs.get("cmap", "Spectral_r")
        im = ax.pcolormesh(x, y, DATA, cmap=cmap, vmin=vmin, vmax=vmax)
        ax.set_aspect("equal")
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.08)
        plt.colorbar(im, ax=ax, cax=cax)
    else:
        raise ValueError("finest_field_plot not yet ready for 3d")

    ax.set_title(kwargs.get("title", ""))

    ax.set_xlabel(kwargs.get("xlabel", "$x \\ (d_p)$"))
    ax.set_ylabel(kwargs.get("ylabel", "$y \\ (d_p)$"))

    if "xlim" in kwargs:
        ax.set_xlim(kwargs["xlim"])
    if "ylim" in kwargs:
        ax.set_ylim(kwargs["ylim"])

    if "filename" in kwargs:
        fig.savefig(kwargs["filename"])

    return fig, ax
