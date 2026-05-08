#
#
#


class PatchLevel:
    """is a collection of patches"""

    def __init__(self, lvl_nbr, patches):
        self.level_number = lvl_nbr
        self.patches = patches

    def __iter__(self):
        return self.patches.__iter__()

    def level_range(self):
        name = list(self.patches[0].patch_datas.keys())[0]
        return min([patch.patch_datas[name].x.min() for patch in self.patches]), max(
            [patch.patch_datas[name].x.max() for patch in self.patches]
        )

    def __getitem__(self, idx):
        if type(idx) is int:
            return self.patches[idx]
        raise IndexError(f"PatchLevel::__getitem__ unhandled input type: {type(idx)}")

    @property
    def cell_width(self):
        return self.patches[0].layout.dl

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        return patch_level_array_ufunc(self, ufunc, method, *inputs, **kwargs)

    def __array_function__(self, func, types, args, kwargs):
        return patch_level_array_function(self, func, types, args, kwargs)


def patch_level_array_ufunc(patch_level, ufunc, method, *inputs, **kwargs):
    if method != "__call__":
        return NotImplemented

    def extract(pidx):
        return [x.patches[pidx] if type(x) is type(patch_level) else x for x in inputs]

    out = [
        getattr(ufunc, method)(*extract(pidx), **kwargs)
        for pidx in range(len(patch_level.patches))
    ]
    return type(patch_level)(patch_level.level_number, out)


def patch_level_array_function(patch_level, func, types, args, kwargs):
    def extract(pidx):
        return [x.patches[pidx] if type(x) is type(patch_level) else x for x in args]

    out = [func(*extract(pidx), **kwargs) for pidx in range(len(patch_level.patches))]

    if type(out[0]) is not type(patch_level[0]):
        return out
    return type(patch_level)(patch_level.level_number, out)
