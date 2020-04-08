def intralevel(diag, patch_level):
    patches = list(patch_level.patches.values())

    if len(patches) == 0:
        return {}

    if diag.dim == 1:
        fn = _border_patches_1d

    if fn is None:
        raise ValueError(
            "Periodic Overlap dimension unsupported, dimension: " + diag.dim
        )

    return {
        data_name: fn(diag, patch_level, data_name)
        for data_name in patch_level.data_names()
    }


def _border_patches_1d(diag, patch_level, data_name):
    patches = list(patch_level.patches.values())
    assert len(patches)

    direction = "x"
    max_domain = diag.domain_upper
    max_x = max_domain[0]
    min_x = diag.origin[0]

    nGhosts = patch_level.nGhosts(data_name)
    ghost_width = patch_level.cell_width(direction) * nGhosts
    ghost_box_lower = min_x + ghost_width
    ghost_box_upper = max_x - ghost_width

    def on_lower_border(patch):
        return patch.min_coord(direction) <= ghost_box_lower

    def on_upper_border(patch):
        return patch.max_coord(direction) >= ghost_box_upper

    def add_if_bordering(on_border, offset=[0, 0, 0]):
        """offset passed to patch.copy to transform patch origin if on_lower_border"""
        return [patch.copy(offset) for patch in patches if on_border(patch)]

    return (
        add_if_bordering(
            on_lower_border,
            [max - min for min, max in zip(diag.origin, diag.domain_upper)],
        ),
        add_if_bordering(on_upper_border),
    )
