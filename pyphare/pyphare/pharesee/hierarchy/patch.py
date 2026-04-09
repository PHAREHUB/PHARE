#

from typing import Any


class Patch:
    """
    A patch represents a hyper-rectangular region of space
    """

    def __init__(self, patch_datas, patch_id="", layout=None, attrs=None):
        """
        :param patch_datas: a list of PatchData objects
        these are assumed to "belong" to the Patch so to
        share the same origin, mesh size and box.
        """
        if layout is not None:
            self.layout = layout
            self.box = layout.box
            self.origin = layout.origin
            self.dl = layout.dl
            self.patch_datas = patch_datas
            self.id = patch_id

        if len(patch_datas):
            pdata0 = list(patch_datas.values())[0]  # 0 represents all others
            self.layout = pdata0.layout
            self.box = pdata0.layout.box
            self.origin = pdata0.layout.origin
            self.dl = pdata0.layout.dl
            self.patch_datas = patch_datas
            self.id = patch_id

        self.attrs = attrs

    def __str__(self):
        return f"Patch: box( {self.box}), id({self.id})"

    def __repr__(self):
        return self.__str__()

    def __iter__(self):
        return self.patch_datas.items().__iter__()

    def __getitem__(self, key):
        if type(key) is Any:
            return next(iter(self.patch_datas.values()))
        if type(key) is int:
            if len(self.patch_datas) != 1:
                raise ValueError(
                    "too many patch_datas to access by int, only allowed if len(patch_datas)== 1",
                    self.patch_datas.keys(),
                )
            return next(iter(self.patch_datas.values()))
        if key in self.patch_datas:
            return self.patch_datas[key]
        for k in self.patch_datas.keys():
            if k.endswith(key):
                return self.patch_datas[k]
        raise ValueError("no key in set", key, self.patch_datas.keys())

    def copy(self):
        """does not copy patchdatas.datasets (see class PatchData)"""
        from copy import deepcopy

        return deepcopy(self)

    def __copy__(self):
        return self.copy()

    def __call__(self, qty, **kwargs):
        # take slice/slab of 1/2d array from 2/3d array
        if "x" in kwargs and len(kwargs) == 1:
            cut = kwargs["x"]
            idim = 0
        elif "y" in kwargs and len(kwargs) == 1:
            cut = kwargs["y"]
            idim = 1
        else:
            raise ValueError("need to specify either x or y cut coordinate")
        pd = self.patch_datas[qty]
        origin = pd.origin[idim]
        idx = int((cut - origin) / pd.layout.dl[idim])
        nbrGhosts = pd.ghosts_nbr[idim]
        if idim == 0:
            return pd.dataset[idx + nbrGhosts, nbrGhosts:-nbrGhosts]
        elif idim == 1:
            return pd.dataset[nbrGhosts:-nbrGhosts, idx + nbrGhosts]

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        return patch_array_ufunc(self, ufunc, method, *inputs, **kwargs)

    def __array_function__(self, func, types, args, kwargs):
        return patch_array_function(self, func, types, args, kwargs)

    def is_compatible(self, that):
        return type(self) is type(that) and self.patch_datas == that.patch_datas


def patch_array_ufunc(patch, ufunc, method, *inputs, **kwargs):
    if method != "__call__":
        return NotImplemented

    if not all([patch.is_compatible(o) for o in inputs]):
        raise TypeError("Patch: incompatible arguments")

    pds = [x.patch_datas.values() for x in inputs]
    out = [getattr(ufunc, method)(*pd, **kwargs) for pd in zip(*pds)]
    final = {k: pd for k, pd in zip(patch.patch_datas, out)}

    return type(patch)(final, patch_id=patch.id, layout=patch.layout, attrs=patch.attrs)


def patch_array_function(patch, func, types, args, kwargs):
    if not all([patch.is_compatible(o) for o in args]):
        raise TypeError("Patch: incompatible arguments")

    pds = [x.patch_datas.values() for x in args]
    out = [func(*pd, **kwargs) for pd in zip(*pds)]
    final = {k: pd for k, pd in zip(patch.patch_datas, out)}
    any_patch_data = next(iter(patch.patch_datas.values()))

    if type(any_patch_data) is not type(out[0]):
        return final

    return type(patch)(final, patch_id=patch.id, layout=patch.layout, attrs=patch.attrs)
