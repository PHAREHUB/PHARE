from .hierarchy import PatchHierarchy


def is_amr_field(obj):
    """True for a field of pyphare.pharesee.amr (detected without importing it): the
    operators of the legacy classes then return NotImplemented and let it answer"""
    return any(
        c.__module__.startswith("pyphare.pharesee.amr") for c in type(obj).__mro__
    )


class AnyTensorField(PatchHierarchy):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @staticmethod
    def FROM(klass, hier):
        return klass(
            hier.time_hier,
            hier.domain_box,
            hier.data_files,
            hier.selection_box,
        )

    def __getitem__(self, input):
        if input in self.__dict__:
            return self.__dict__[input]
        raise IndexError("AnyTensorField.__getitem__ cannot handle input", input)
