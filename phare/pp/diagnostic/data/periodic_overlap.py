def intralevel(diag, patch_level):
    patches = list(patch_level.patches.values())

    if len(patches) == 0:
        return {}

    if diag.dim == 1:
        border_patches = _border_patches_1d

    if border_patches is None:
        raise ValueError(
            "Periodic Overlap dimension unsupported, dimension: " + diag.dim
        )

    return {
        data_name: border_patches(diag, patch_level, data_name)
        for data_name in patch_level.data_names()
    }


def _border_patches_1d(diag, patch_level, data_name):
    patches = sorted(patch_level.patches.values(), key=lambda x: x.min_coord("x"))
    assert len(patches)

    direction = "x"
    min_x, max_x = diag.origin[0], diag.domain_upper[0]
    ghost_width = patch_level.cell_width(direction) * patch_level.nGhosts(data_name)
    ghost_box_lower, ghost_box_upper = min_x + ghost_width, max_x - ghost_width
    on_lower_border = patches[0].min_coord(direction) <= ghost_box_lower
    on_upper_border = patches[-1].max_coord(direction) >= ghost_box_upper
    lower_offset = [max - min for min, max in zip(diag.origin, diag.domain_upper)]

    return (
        [patches[0].copy(lower_offset)] if on_lower_border else [],
        [patches[-1]] if on_upper_border else [],
    )
