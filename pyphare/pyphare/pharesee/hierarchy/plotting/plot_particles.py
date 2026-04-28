def hierarchy_dist_plot(hier, **kwargs):
    import copy
    import numpy as np
    from pyphare.core import phare_utilities as phut

    usr_lvls = kwargs.get("levels", (0,))
    finest = kwargs.get("finest", False)
    pops = kwargs.get("pop", [])
    time = kwargs.get("time", hier.times()[0])
    axis = phut.listify(kwargs.get("axis", ("Vx", "Vy")))
    all_pops = list(hier.level(0, time).patches[0].patch_datas.keys())

    vmin = kwargs.get("vmin", -2)
    vmax = kwargs.get("vmax", 2)
    dv = kwargs.get("dv", 0.05)
    vbins = vmin + dv * np.arange(int((vmax - vmin) / dv))

    if finest:
        bins = None
        final = finest_part_data(hier, time)
        if len(axis) == 2:
            if axis[0] == "x":
                xbins = amr_grid(hier, time)
                bins = (xbins, vbins)
            else:
                bins = (vbins, vbins)
        elif len(axis) == 1:
            bins = vbins
        kwargs["bins"] = bins

    else:
        final = {pop: None for pop in all_pops}
        for lvl_nbr, level in hier.levels(time).items():
            if lvl_nbr not in usr_lvls:
                continue
            for ip, patch in enumerate(level.patches):
                if len(pops) == 0:
                    pops = list(patch.patch_datas.keys())

                for pop in pops:
                    if patch.patch_datas.__len__() == 0:
                        tmp = None
                    else:
                        tmp = copy.copy(patch.patch_datas[pop].dataset)

                    if tmp is not None:
                        if final[pop] is None:
                            final[pop] = tmp
                        else:
                            final[pop].add(tmp)

    # select particles
    if "select" in kwargs:
        for pop, particles in final.items():
            final[pop] = kwargs["select"](particles)

    return final, particles_dist_plot(final, **kwargs)


def finest_part_data(hierarchy, time=None):
    """
    returns a dict {popname : Particles}
    Particles contained in the dict are those from
    the finest patches available at a given location
    """
    import numpy as np
    from copy import deepcopy
    from pyphare.core import box as boxm
    from pyphare.core import phare_utilities as phut
    from pyphare.pharesee.particles import remove

    i_ref = None
    for i_ref, p in enumerate(hierarchy.level(0, time=time).patches):
        if len(p.patch_datas):
            break
    if i_ref is None:
        raise ValueError("This particle hierarchy seems empty !")

    # we are going to return a dict {popname : Particles}
    # we prepare it with population names
    aPatch = hierarchy.level(0, time=time).patches[i_ref]
    particles = {popname: None for popname in aPatch.patch_datas.keys()}

    # our strategy is to explore the hierarchy from the finest
    # level to the coarsest. at Each level we keep only particles
    # that are in cells that are not overlaped by finer boxes

    # this dict keeps boxes for patches at each level
    # each level will thus need this dict to see next finer boxes
    lvlPatchBoxes = {ilvl: [] for ilvl in range(hierarchy.finest_level(time) + 1)}

    for ilvl in range(hierarchy.finest_level(time) + 1)[::-1]:
        plvl = hierarchy.level(ilvl, time=time)
        for ip, patch in enumerate(plvl.patches):
            lvlPatchBoxes[ilvl].append(patch.box)
            for popname, pdata in patch.patch_datas.items():
                # if we're at the finest level
                # we need to keep all particles
                if ilvl == hierarchy.finest_level(time):
                    if particles[popname] is None:
                        particles[popname] = deepcopy(pdata.dataset)
                    else:
                        particles[popname].add(deepcopy(pdata.dataset))

                # if there is a finer level
                # we need to keep only those of the current patch
                # that are not in cells overlaped by finer patch boxes
                else:
                    icells = pdata.dataset.iCells
                    parts = deepcopy(pdata.dataset)
                    create = True
                    for finerBox in lvlPatchBoxes[ilvl + 1]:
                        coarseFinerBox = boxm.coarsen(finerBox, phut.refinement_ratio)
                        within = np.where(
                            (icells >= coarseFinerBox.lower[0])
                            & (icells <= coarseFinerBox.upper[0])
                        )[0]
                        if create:
                            toRemove = within
                            create = False
                        else:
                            toRemove = np.concatenate((toRemove, within))

                    if toRemove.size != 0:
                        parts = remove(parts, toRemove)
                    if parts is not None:
                        particles[popname].add(parts)
    return particles


def amr_grid(hierarchy, time):
    """returns a non-uniform contiguous primal grid
    associated to the given hierarchy
    """
    import numpy as np
    from pyphare.core import box as boxm
    from pyphare.core import phare_utilities as phut

    lvlPatchBoxes = {ilvl: [] for ilvl in range(hierarchy.finest_level() + 1)}
    finalCells = {ilvl: None for ilvl in range(hierarchy.finest_level() + 1)}
    lvl = hierarchy.levels(time)

    for ilvl in range(hierarchy.finest_level(time) + 1)[::-1]:
        sorted_patches = sorted(lvl[ilvl].patches, key=lambda p: p.layout.box.lower[0])

        for ip, patch in enumerate(sorted_patches):
            box = patch.layout.box
            lvlPatchBoxes[ilvl].append(box)

            # we create a list of all cells in the current patch
            # remember that if the box upper cell is, say = 40,
            # it means that the upper node is the lower node of cell 41
            # so to get all primal nodes of a patch we need to include
            # one past the upper cell.
            # this said we do not want to include that last primal nodes
            # all the time because that would be a duplicate with the lower
            # node of the next patch. We only want to add it for the LAST
            # (because sorted) patch. We also do not want to do it on levels
            # other than L0 because the last primal node of the last patch
            # of L_i is the first primal node of a L_{i-1} node, so including it
            # would also mean adding a duplicate.
            last = 1 if ilvl == 0 and ip == len(sorted_patches) - 1 else 0
            cells = np.arange(box.lower[0], box.upper[0] + 1 + last)

            # finest level has no next finer so we take all cells
            if ilvl == hierarchy.finest_level(time):
                if finalCells[ilvl] is None:
                    finalCells[ilvl] = cells
                else:
                    finalCells[ilvl] = np.concatenate((finalCells[ilvl], cells))

            else:
                # on other levels
                # we take only grids not overlaped by next finer
                coarsenedNextFinerBoxes = [
                    boxm.coarsen(b, phut.refinement_ratio)
                    for b in lvlPatchBoxes[ilvl + 1]
                ]
                for coarseBox in coarsenedNextFinerBoxes:
                    ccells = np.arange(coarseBox.lower[0], coarseBox.upper[0] + 1)
                    inter, icells, iccells = np.intersect1d(
                        cells, ccells, return_indices=True
                    )
                    cells = np.delete(cells, icells)
                if len(cells):
                    if finalCells[ilvl] is None:
                        finalCells[ilvl] = cells
                    else:
                        finalCells[ilvl] = np.unique(
                            np.concatenate((finalCells[ilvl], cells))
                        )

    # now we have all cells for each level we
    # just need to compute the primal coordinates
    # and concatenate in a single array
    for ilvl in range(hierarchy.finest_level() + 1):
        if ilvl == 0:
            x = finalCells[ilvl] * hierarchy.level(ilvl).patches[0].layout.dl[0]
        else:
            xx = finalCells[ilvl] * hierarchy.level(ilvl).patches[0].layout.dl[0]
            x = np.concatenate((x, xx))

    return np.sort(x)


def particles_dist_plot(particles, **kwargs):
    """
    plot the phase space of given particles
    particles can be of type Particles, list(Particles), dict{popname:Particles}

    kwargs:
    * axis : ("x", "Vx"), ("x", "Vy"), ("x", "Vz"), ("Vx", "Vy") (default) --
       ("Vx", "Vz"), ("Vy", "Vz"), ("Vx"), ("Vy"), ("Vz")
    * bins :  number of bins in each dimension, default is (50,50)
    * gaussian_filter_sigma : sigma of the gaussian filter, default is (0,0)
    * median_filter_size : size of the median filter, default is (0,0)
    * cmap : color table, default is "jet"
    * norm  : histogram will be normed to Normalize(0,norm)
    * kde : (default False) : adds contours of kernel density estimate
    * title : (str) title of the plot
    * xlabel, ylabel
    * xlim, ylim
    * bulk : (bool) (default : False), adds vertical/horizontal lines --
             at in-plane bulk velocity for velocity axis
    * filename : (str) if exists, save plot to figure under that name

    return value : fig,ax
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm, Normalize
    from pyphare.pharesee.particles import Particles, aggregate
    from pyphare.core.phare_utilities import listify
    from scipy.interpolate import LinearNDInterpolator

    if isinstance(particles, list):
        particles = aggregate(particles)
    elif isinstance(particles, dict):
        particles = aggregate([p for p in particles.values()])

    if not isinstance(particles, Particles):
        raise ValueError("Error, 'particles' type should be Particles, list or dict")

    if "ax" not in kwargs:
        fig, ax = plt.subplots()
    else:
        ax = kwargs["ax"]
        fig = ax.figure
    axis = listify(kwargs.get("axis", ("Vx", "Vy")))
    vaxis = {"Vx": 0, "Vy": 1, "Vz": 2}

    if len(axis) == 2:
        if axis[0] in vaxis:
            x = particles.v[:, vaxis[axis[0]]]
        elif axis[0] == "x":
            x = particles.x
        else:
            raise ValueError("Only abscissa and velocity X-axis are supported yet")
        if axis[1] in vaxis:
            y = particles.v[:, vaxis[axis[1]]]
        else:
            raise ValueError("Only velocity Y-axis are supported yet")

        bins = kwargs.get("bins", (50, 50))
        h, xh, yh = np.histogram2d(x, y, bins=bins, weights=particles.weights[:, 0])

        if "gaussian_filter_sigma" in kwargs and "median_filter_size" not in kwargs:
            from scipy.ndimage import gaussian_filter

            sig = kwargs.get("gaussian_filter_sigma", (0, 0))
            image = gaussian_filter(h.T, sigma=sig)
        elif "median_filter_size" in kwargs and "gaussian_filter_sigma" not in kwargs:
            from scipy.ndimage import median_filter

            siz = kwargs.get("median_filter_size", (0, 0))
            image = median_filter(h.T, size=siz)
        elif (
            "gaussian_filter_sigma" not in kwargs and "median_filter_size" not in kwargs
        ):
            image = h.T
        else:
            raise ValueError(
                "gaussian and median filters can not be called at the same time"
            )

        plain = kwargs.get("plain", False)

        if not plain:
            cmap = kwargs.get("cmap", "jet")

            cmax = kwargs.get("color_max", h.max())
            cmin = kwargs.get("color_min", h.min())
            cmin = max(cmin, 1e-4)

            color_scale = kwargs.get("color_scale", "log")
            if color_scale == "log":
                norm = LogNorm(vmin=cmin, vmax=cmax)
            elif color_scale == "linear":
                norm = Normalize(cmin, cmax)
            else:
                raise ValueError(
                    "Only log and linear color_scale values are supported yet"
                )

            im = ax.pcolormesh(xh, yh, image, cmap=cmap, norm=norm)

            fig.colorbar(im, ax=ax)
        else:
            color = kwargs.get("color", "k")
            stride = kwargs.get("stride", 1)
            if stride <= 0:
                raise ValueError("stride must be a positive integer")
            markersize = kwargs.get("markersize", 0.5)
            alpha = kwargs.get("alpha", 0.5)
            im = ax.scatter(
                x[::stride],
                y[::stride],
                color=color,
                marker=".",
                markersize=markersize,
                alpha=alpha,
            )

        ax.set_ylabel(kwargs.get("ylabel", axis[1]))

    elif len(axis) == 1:
        bins = kwargs.get("bins", (50))
        cuts = kwargs.get("cuts", None)
        ndim = particles.ndim

        if cuts is not None:
            from pyphare.core.box import Box

            if ndim == 1:
                if cuts[0] is None:
                    new_particles = particles
                else:
                    box_new = Box(cuts[0][0], cuts[0][1])
                    new_particles = particles.select(box_new, box_type="pos")
            elif ndim == 2:
                box_new = Box(
                    (cuts[0][0], cuts[1][0]), (cuts[0][1], cuts[1][1])
                )  # TODO need to be tested
                new_particles = particles.select(box_new, box_type="pos")
            else:
                box_new = Box(
                    (cuts[0][0], cuts[1][0], cuts[2][0]),
                    (cuts[0][1], cuts[1][1], cuts[2][1]),
                )  # TODO need to be tested
                new_particles = particles.select(box_new, box_type="pos")
        else:
            new_particles = particles

        drawstyle = kwargs.get("drawstyle", "steps-mid")

        if axis[0] in vaxis:
            x = new_particles.v[:, vaxis[axis[0]]]
        else:
            raise ValueError(
                "For 1-D dist_plot, the abscissa has to be a velocity axis"
            )

        bins = kwargs.get("bins", 50)
        h, xh = np.histogram(x, bins=bins, weights=new_particles.weights[:, 0])

        xh_ = 0.5 * (xh[:-1] + xh[1:])
        im = ax.plot(xh_, h, drawstyle=drawstyle)

        ax.set_ylabel(kwargs.get("ylabel", "f"))

    if kwargs.get("kde", False) is True:
        import seaborn as sns

        if len(axis) == 1:
            sns.kdeplot(x=x, ax=ax, color="w")
        else:
            sns.kdeplot(x=x, y=y, ax=ax, color="w")

    ax.set_title(kwargs.get("title", ""))
    ax.set_xlabel(kwargs.get("xlabel", axis[0]))

    if "xlim" in kwargs:
        ax.set_xlim(kwargs["xlim"])
    if "ylim" in kwargs:
        ax.set_ylim(kwargs["ylim"])

    if "bulk" in kwargs:
        if kwargs["bulk"] is True:
            if axis[0] in vaxis:
                ax.axvline(
                    np.average(
                        new_particles.v[:, vaxis[axis[0]]],
                        weights=new_particles.weights,
                    ),
                    color="w",
                    ls="--",
                )
            if axis[1] in vaxis:
                ax.axhline(
                    np.average(
                        new_particles.v[:, vaxis[axis[1]]],
                        weights=new_particles.weights,
                    ),
                    color="w",
                    ls="--",
                )

    if "filename" in kwargs:
        fig.savefig(kwargs["filename"])

    interp = kwargs.get("interp", False)

    if interp:
        xbins = 0.5 * (xh[1:] + xh[:-1])
        ybins = 0.5 * (yh[1:] + yh[:-1])
        xx, yy = np.meshgrid(xbins, ybins, indexing="ij")
        coords = np.array([xx.flatten(), yy.flatten()]).T
        interpdist = LinearNDInterpolator(coords, image.T.flatten())
        return fig, ax, interpdist, xbins, ybins

    return fig, ax


def connect_bbox(
    bbox1, bbox2, loc1a, loc2a, loc1b, loc2b, prop_lines, prop_patches=None
):
    from mpl_toolkits.axes_grid1.inset_locator import (
        BboxPatch,
        BboxConnector,
        BboxConnectorPatch,
    )

    if prop_patches is None:
        prop_patches = {
            **prop_lines,
            "alpha": prop_lines.get("alpha", 1) * 0.2,
        }

    c1 = BboxConnector(bbox1, bbox2, loc1=loc1a, loc2=loc2a, **prop_lines)
    c1.set_clip_on(False)
    c2 = BboxConnector(bbox1, bbox2, loc1=loc1b, loc2=loc2b, **prop_lines)
    c2.set_clip_on(False)
    bbox_patch1 = BboxPatch(bbox1, ec="k", fc="none", ls="--")
    bbox_patch2 = BboxPatch(bbox2, ec="k", fc="none", ls="--")

    p = BboxConnectorPatch(
        bbox1,
        bbox2,
        # loc1a=3, loc2a=2, loc1b=4, loc2b=1,
        loc1a=loc1a,
        loc2a=loc2a,
        loc1b=loc1b,
        loc2b=loc2b,
        **prop_patches
    )
    p.set_clip_on(False)

    return c1, c2, bbox_patch1, bbox_patch2, p


def zoom_effect(ax1, ax2, xmin, xmax, **kwargs):
    """
    Connect *ax1* and *ax2*. The *xmin*-to-*xmax* range in both axes will
    be marked.

    Parameters
    ----------
    ax1
        The main axes.
    ax2
        The zoomed axes.
    xmin, xmax
        The limits of the colored area in both plot axes.
    **kwargs
        Arguments passed to the patch constructor.
    """
    from matplotlib.transforms import Bbox, TransformedBbox, blended_transform_factory

    trans1 = blended_transform_factory(ax1.transData, ax1.transAxes)
    trans2 = blended_transform_factory(ax2.transData, ax2.transAxes)

    bbox = Bbox.from_extents(xmin, 0, xmax, 1)

    mybbox1 = TransformedBbox(bbox, trans1)
    mybbox2 = TransformedBbox(bbox, trans2)

    prop_patches = {**kwargs, "ec": "none", "alpha": 0.2}

    c1, c2, bbox_patch1, bbox_patch2, p = connect_bbox(
        mybbox1,
        mybbox2,
        loc1a=3,
        loc2a=2,
        loc1b=4,
        loc2b=1,
        prop_lines=kwargs,
        prop_patches=prop_patches,
    )

    ax1.add_patch(bbox_patch1)
    ax2.add_patch(bbox_patch2)
    ax2.add_patch(c1)
    ax2.add_patch(c2)
    ax2.add_patch(p)

    return c1, c2, bbox_patch1, bbox_patch2, p
