#
#
#


from . import tensorfield
from .hierarchy_utils import compute_hier_from, compute_rename, rename, _compute_neg


class ScalarField(tensorfield.AnyTensorField):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @classmethod
    def FROM(cls, hier):
        renamed_hier = compute_hier_from(compute_rename, hier, new_names=("value",))
        return super().FROM(cls, renamed_hier)

    def __add__(self, other):
        assert isinstance(other, (ScalarField, int, float))
        h_self = rename(self, ["self_value"])

        if isinstance(other, ScalarField):
            h_other = rename(other, ["other_value"])
            h = compute_hier_from(
                self._compute_add,
                (h_self, h_other),
            )
        elif isinstance(other, (int, float)):
            h = compute_hier_from(self._compute_add, (h_self,), other=other)
        else:
            raise RuntimeError("right operand not supported")

        return ScalarField.FROM(h)

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        assert isinstance(other, (ScalarField, int, float))
        h_self = rename(self, ["self_value"])

        if isinstance(other, ScalarField):
            h_other = rename(other, ["other_value"])
            h = compute_hier_from(
                self._compute_sub,
                (h_self, h_other),
            )
        elif isinstance(other, (int, float)):
            h = compute_hier_from(self._compute_sub, (h_self,), other=other)
        else:
            raise RuntimeError("right operand not supported")

        return ScalarField.FROM(h)

    def __mul__(self, other):
        assert isinstance(other, (ScalarField, int, float))
        h_self = rename(self, ["self_value"])

        if isinstance(other, ScalarField):
            h_other = rename(other, ["other_value"])
            h = compute_hier_from(self._compute_mul, (h_self, h_other))
        elif isinstance(other, (int, float)):
            h = compute_hier_from(self._compute_mul, (h_self,), other=other)
        else:
            raise RuntimeError("right operand not supported")

        return ScalarField.FROM(h)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        assert isinstance(other, (ScalarField, int, float))
        h_self = rename(self, ["self_value"])

        if isinstance(other, ScalarField):
            h_other = rename(other, ["other_value"])
            h = compute_hier_from(self._compute_truediv, (h_self, h_other))
        elif isinstance(other, (int, float)):
            h = compute_hier_from(self._compute_truediv, (h_self,), other=other)
        else:
            raise RuntimeError("right operand not supported")

        return ScalarField.FROM(h)

    def __rtruediv__(self, other):
        assert isinstance(other, (int, float))
        h_self = rename(self, ["self_value"])

        h = compute_hier_from(self._compute_rtruediv, (h_self,), other=other)

        return ScalarField.FROM(h)

    def _compute_add(self, patch_datas, **kwargs):
        ref_name = next(iter(patch_datas.keys()))

        if "other" in kwargs:
            other = kwargs["other"]
            dset_value = patch_datas["self_value"].dataset[:] + other
        else:
            dset_value = (
                patch_datas["self_value"].dataset[:]
                + patch_datas["other_value"].dataset[:]
            )

        return (
            {
                "name": "value",
                "data": dset_value,
                "centering": patch_datas[ref_name].centerings,
            },
        )

    def _compute_sub(self, patch_datas, **kwargs):
        # subtract a scalar from the dataset of a scalarField
        if "other" in kwargs:
            other = kwargs["other"]
            dset_value = patch_datas["self_value"].dataset[:] - other
        # substraction of 2 scalarFields
        else:
            dset_value = (
                patch_datas["self_value"].dataset[:]
                - patch_datas["other_value"].dataset[:]
            )

        return (
            {
                "name": "value",
                "data": dset_value,
                "centering": patch_datas["self_value"].centerings,
            },
        )

    def _compute_mul(self, patch_datas, **kwargs):
        # multiplication of a scalarField by a scalar
        if "other" in kwargs:
            other = kwargs["other"]
            pd_attrs = []

            for pd_name in patch_datas:
                pd_attrs.append(
                    {
                        "name": "value",
                        "data": other * patch_datas[pd_name].dataset[:],
                        "centering": patch_datas[pd_name].centerings,
                    }
                )
        # multiplication of 2 scalarField
        else:
            dset_value = (
                patch_datas["self_value"].dataset[:]
                * patch_datas["other_value"].dataset[:]
            )
            pd_attrs = (
                {
                    "name": "value",
                    "data": dset_value,
                    "centering": patch_datas["self_value"].centerings,
                },
            )

        return tuple(pd_attrs)

    def _compute_truediv(self, patch_datas, **kwargs):
        # multiplication of a scalarField by a scalar
        if "other" in kwargs:
            other = kwargs["other"]
            pd_attrs = []

            for pd_name in patch_datas:
                pd_attrs.append(
                    {
                        "name": "value",
                        "data": patch_datas[pd_name].dataset[:] / other,
                        "centering": patch_datas[pd_name].centerings,
                    }
                )
        # multiplication of 2 scalarField
        else:
            dset_value = (
                patch_datas["self_value"].dataset[:]
                / patch_datas["other_value"].dataset[:]
            )
            pd_attrs = (
                {
                    "name": "value",
                    "data": dset_value,
                    "centering": patch_datas["self_value"].centerings,
                },
            )

        return tuple(pd_attrs)

    def _compute_rtruediv(self, patch_datas, **kwargs):
        # Scalar divided by a scalarField
        other = kwargs["other"]
        pd_attrs = []

        for pd_name in patch_datas:
            pd_attrs.append(
                {
                    "name": "value",
                    "data": other / patch_datas[pd_name].dataset[:],
                    "centering": patch_datas[pd_name].centerings,
                }
            )

        return tuple(pd_attrs)

    def __neg__(self):
        names_self = self.quantities()
        h = compute_hier_from(_compute_neg, self, new_names=names_self)
        return ScalarField.FROM(h)
