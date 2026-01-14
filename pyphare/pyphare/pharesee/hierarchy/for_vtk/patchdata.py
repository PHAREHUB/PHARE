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


class VtkFieldData(patchdata.FieldData):
    def __init__(self, lvl_info, patch_idx, cmp_idx, data_offset, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.lvl_info = lvl_info
        self.patch_idx = patch_idx
        self.cmp_idx = cmp_idx
        self.data_offset = data_offset
        self.dataset = (
            self.dataset
            if type(self.dataset) is VtkFieldDatasetAccessor
            else VtkFieldDatasetAccessor(self, self.dataset)
        )
        self.field_box = boxm.Box(self.box.lower, self.box.upper + 1)
        self.data_size = self.field_box.size()

    def compare(self, that, atol=1e-16):
        """VTK Diagnostics do not have ghosts values!"""
        try:
            that_data = that[:] if type(that) is VtkFieldData else that[that.box]
            phut.assert_fp_any_all_close(self.dataset[:], that_data, atol=atol)
            return True
        except Exception as e:
            return phut.EqualityCheck(False, str(e))

    def __eq__(self, that):
        return self.compare(that)

    def __deepcopy__(self, memo):
        no_copy_keys = ["dataset"]  # do not copy these things
        cpy = phut.deep_copy(self, memo, no_copy_keys)
        cpy.dataset = self.dataset
        return cpy
