#
#
#


from pyphare.pharesee.hierarchy import patch


class Patch(patch.Patch):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def __getitem__(self, key):
        return self.patch_datas[key]

    def copy(self):
        from copy import deepcopy

        return deepcopy(self)

    def __copy__(self):
        return self.copy()

    def __call__(self, qty, **kwargs):
        raise RuntimeError("finish")
