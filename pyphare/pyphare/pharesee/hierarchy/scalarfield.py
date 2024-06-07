from .hierarchy import PatchHierarchy
from .hierarchy_utils import compute_hier_from, compute_rename, rename, _compute_neg


class ScalarField(PatchHierarchy):
    def __init__(self, hier):
        renamed_hier = compute_hier_from(compute_rename, hier, new_names=("value",))
        patch_levels = renamed_hier.patch_levels
        domain_box = renamed_hier.domain_box
        refinement_ratio = renamed_hier.refinement_ratio
        data_files = renamed_hier.data_files

        super().__init__(
            patch_levels, domain_box, refinement_ratio, renamed_hier.times(), data_files
        )

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

        return ScalarField(h)

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

        return ScalarField(h)

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

        return ScalarField(h)

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

        return ScalarField(h)

    def __rtruediv__(self, other):
        assert isinstance(other, (int, float))
        h_self = rename(self, ["self_value"])

        h = compute_hier_from(self._compute_rtruediv, (h_self,), other=other)

        return ScalarField(h)

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
        return ScalarField(h)
