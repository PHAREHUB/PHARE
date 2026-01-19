#
#
#


import pyphare.core.box as boxm
from pyphare.pharesee.hierarchy import patchdata

from pyphare.core import phare_utilities as phut


class VtkFieldDatasetAccessor:
    def __init__(self, field_data, dataset):
        self.field_data = field_data
        self.dataset = dataset

    def __getitem__(self, slice):  # todo finish slice/box lookup
        cmp_idx = self.field_data.cmp_idx
        off_set = self.field_data.data_offset
        data_size = self.field_data.data_size
        return self.dataset[:, cmp_idx][off_set : off_set + data_size].reshape(
            self.field_data.field_box.shape, order="F"
        )

    @property
    def shape(self):
        return self.field_data.field_box.shape


class VtkFieldData(patchdata.FieldData):
    def __init__(self, layout, cmp_idx, data_offset, *args, **kwargs):
        super().__init__(
            layout,
            *args,
            centering=["primal"] * layout.box.ndim,
            ghosts_nbr=[0] * layout.box.ndim,
            **kwargs
        )

        import h5py

        self.cmp_idx = cmp_idx
        self.data_offset = data_offset
        self.dataset = (
            VtkFieldDatasetAccessor(self, self.dataset)
            if type(self.dataset) is h5py.Dataset
            else self.dataset
        )
        self.field_box = boxm.Box(self.box.lower, self.box.upper + 1)
        self.data_size = self.field_box.size()

    def compare(self, that, atol=1e-16):
        """VTK Diagnostics do not have ghosts values!"""

        try:
            that_data = (
                that[:]
                if all([that.dataset.shape == self.dataset.shape])
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
        return type(self)(
            self.layout,
            self.cmp_idx,
            self.data_offset,
            self.field_name,
            data if data is not None else self.dataset,
            **kwargs
        )
