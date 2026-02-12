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

        # Collect per-argument sequences
        seqs = []
        for x in inputs:
            print("zobi : ", type(x))
            if isinstance(x, Patch):
                print("zob", x.patch_datas)
                for k, v in x.patch_datas.items():
                    print("toto", k, v, type(v))
                    seqs.append(v)
                    if isinstance(v, FieldData):
                        print("is a fieldData")
                #seqs.append(x.patch_datas)
            else:
                #  TODO seqs.append(x) ?
                raise TypeError("this arg should be a Patch")

        # print(type(seqs[0]))
        # print(seqs[0].dataset.shape)
        # # Outer length must match
        # n = len(seqs[0])                        # length of the first patch
        # if not all(len(s) == n for s in seqs):  # each of the other Patches have to be homogeneous
        #     raise ValueError("Patch length mismatch")

        print(seqs)
        for elem in zip(*seqs):
            print(elem)

        # Pairwise application
        result = [
            ufunc(*elems, **kwargs)
            for elems in zip(*seqs)
        ]

        return Patch(result, patch_id=self.id, layout=self.layout, attrs=self.attrs)

    def __array_function__(self, func, types, args, kwargs):
        print(f"__array_function__ of Patch {func.__name__} called for {[getattr(a, 'name', a) for a in args]}")
        return func(*args, **kwargs)

