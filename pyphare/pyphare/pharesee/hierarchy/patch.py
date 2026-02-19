#
from .patchdata import FieldData


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

    def __getitem__(self, key):
        return self.patch_datas[key]

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
        print(f"__array_function__ of Patch called for {ufunc.__name__}")
        if method != "__call__":
            return NotImplemented

        pds = []  # list (p1, p2, p3... ) of list ('x', 'y', 'z') of pds
        for x in inputs:  # inputs is a list of Patch
            if isinstance(x, Patch):  # hence, x is a Patch
                pd_k = []
                pds_ = []
                for k, p in x.patch_datas.items():
                    pd_k.append(k)
                    pds_.append(p)
                pds.append(pds_)
            else:
                raise TypeError("this arg should be a Patch")

        out = [getattr(ufunc, method)(*pd, **kwargs)  for pd in zip(*pds)]
        # out = [ufunc(*pd, **kwargs) for pd in zip(*pds)]

        final = {}
        for k, pd in zip(pd_k, out):  # TODO hmmmm, the output patch will keep the keys of the last patch in inputs
            final[k] = pd

        return Patch(final, patch_id=self.id, layout=self.layout, attrs=self.attrs)

    def __array_function__(self, func, types, args, kwargs):
        # TODO this has to be tested w. np.mean for example
        print(f"__array_function__ of Patch {func.__name__} called for {[getattr(a, 'name', a) for a in args]}")
        return func(*args, **kwargs)

