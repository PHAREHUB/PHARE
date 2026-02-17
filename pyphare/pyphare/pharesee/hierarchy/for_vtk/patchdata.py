#
#
#


import pyphare.core.box as boxm
from pyphare.pharesee.hierarchy import patchdata

from pyphare.core import phare_utilities as phut


class VtkFieldDatasetAccessor:
    def __init__(self, dataset, cmp_idx, offset, box):
        self.box = box
        self.cmp_idx = cmp_idx
        self.offset = offset
        self.dataset = dataset

    def __getitem__(self, slice):
        # todo finish slice/box lookup
        return self.dataset[:, self.cmp_idx][
            self.offset : self.offset + self.box.size()
        ].reshape(self.box.shape, order="F")

    @property
    def shape(self):
        return self.box.shape


class VtkFieldData(patchdata.FieldData):
    def __init__(self, layout, name, data, cmp_idx, offset, **kwargs):
        import h5py

        data_t = type(data)
        if not (data_t is h5py.Dataset or data_t is VtkFieldDatasetAccessor):
            raise RuntimeError("VtkFieldData only handles vtkhdf datasets")

        box = layout.box
        self.field_box = boxm.Box(box.lower, box.upper + 1)
        kwargs["ghosts_nbr"] = [0] * layout.box.ndim
        kwargs["centering"] = ["primal"] * layout.box.ndim
        data = (
            data
            if data_t is VtkFieldDatasetAccessor
            else VtkFieldDatasetAccessor(data, cmp_idx, offset, self.field_box)
        )
        super().__init__(layout, name, data, **kwargs)

    def compare(self, that, atol=1e-16):
        """VTK Diagnostics do not have ghosts values!"""

        try:
            that_data = (
                that[:]
                if (that.dataset.shape == self.dataset.shape).all()
                else that[that.box]
            )
            phut.assert_fp_any_all_close(self.dataset[:], that_data, atol=atol)
            return True
        except AssertionError as e:
            return phut.EqualityCheck(False, str(e))

    def __eq__(self, that):
        return self.compare(that)

    def __deepcopy__(self, memo):
        no_copy_keys = ["dataset"]  # do not copy these things
        cpy = phut.deep_copy(self, memo, no_copy_keys)
        cpy.dataset = self.dataset
        return cpy

    def copy_as(self, data=None, **kwargs):
        data = self.dataset if data is None else data
        if type(data) is VtkFieldDatasetAccessor:
            return type(self)(
                self.layout, self.field_name, data, data.cmp_idx, data.offset, **kwargs
            )
        # make a normal FieldData
        return super().copy_as(data, **kwargs)
