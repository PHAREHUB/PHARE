def intralevel(clazz, diag, patches):
    if len(patches) == 0:
        return {}
    if diag.dim == 1:
        fn = _border_patches_1d

    if fn is None:
        raise ValueError(
            "Periodic Overlap dimension unsupported, dimension: " + diag.dim
        )
    return {
        dataset: fn(diag, patches, patches[0].dtype.nGhosts(dataset))
        for dataset in patches[0].dtype.keys()
    }


def _border_patches_1d(diag, patches, nGhosts):
    assert len(patches)

    direction = "x"
    max_domain = diag.sim.simulation_domain()
    max_x = max_domain[0]
    min_x = diag.sim.origin[0]

    ghost_width = patches[0].patch_level.cell_width(direction) * nGhosts
    lower_ghost_width = min_x + ghost_width
    upper_ghost_width = max_x - ghost_width

    def on_lower_border(patch):
        return patch.min_coord(direction) <= lower_ghost_width

    def on_upper_border(patch):
        return patch.max_coord(direction) >= upper_ghost_width

    def add_if_bordering(on_border, offset=[0, 0, 0]):
        """offset passed to patch.copy to transform patch origin if on_lower_border"""
        return [patch.copy(offset) for patch in patches if on_border(patch)]

    return (
        add_if_bordering(on_lower_border, max_domain),
        add_if_bordering(on_upper_border),
    )
