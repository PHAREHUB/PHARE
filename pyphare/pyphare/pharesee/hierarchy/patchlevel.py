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

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        print(f"__array_function__ of PatchLevel called for {ufunc.__name__}")
        if method != "__call__":
            return NotImplemented

        ps = []
        for x in inputs:
            if isinstance(x, PatchLevel):
                ps_ = []
                for p in x.patches:
                    ps_.append(p)
                ps.append(ps_)
            else:
                raise TypeError("this arg should be a PatchLevel")

        out = [getattr(ufunc, method)(*p, **kwargs)  for p in zip(*ps)]

        return PatchLevel(self.level_number, out)

    def __array_function__(self, func, types, args, kwargs):
        # TODO this has to be tested w. np.mean for example
        print(f"__array_function__ of Patch {func.__name__} called for {[getattr(a, 'name', a) for a in args]}")
        return func(*args, **kwargs)

