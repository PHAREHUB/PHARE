from pyphare.pharesee.hierarchy.hierarchy_utils import compute_hier_from


def hierarchy_from_func1d(func, hier, **kwargs):
    assert hier.ndim == 1

    def compute_(patch_datas, **kwargs):
        ref_name = next(iter(patch_datas.keys()))
        x_ = patch_datas[ref_name].x

        return (
            {"name": "value", "data": func(x_, **kwargs), "centering": patch_datas[ref_name].centerings},
        )

    return compute_hier_from(compute_, hier, **kwargs)


def hierarchy_from_func2d(func, hier, **kwargs):
    from pyphare.pharesee.hierarchy.hierarchy_utils import compute_hier_from

    def compute_(patch_datas, **kwargs):
        ref_name = next(iter(patch_datas.keys()))
        x_ = patch_datas[ref_name].x
        y_ = patch_datas[ref_name].y

        return (
            {"name": "value", "data": func(x_, y_, **kwargs), "centering": patch_datas[ref_name].centerings},
        )

    return compute_hier_from(compute_, hier, **kwargs)


def hierarchy_from_func(func, hier, **kwargs):

    """
    Route hierarchical computation to appropriate dimension handler.

    Parameters
    ----------
    func : callable
        Function to apply to coordinates of the hierarchy
    hier : Hierarchy
        Hierarchy object (1D or 2D)
    **kwargs : dict
        Additional arguments passed to func

    Returns
    -------
    dict
        Computed hierarchical data

    Raises
    ------
    ValueError
        If hierarchy dimension is not supported
    """

    if hier.ndim == 1:
        return hierarchy_from_func1d(func, hier, **kwargs)
    if hier.ndim == 2:
        return hierarchy_from_func2d(func, hier, **kwargs)
    else:
        raise ValueError(f"Unsupported hierarchy dimension: {hier.ndim}")
