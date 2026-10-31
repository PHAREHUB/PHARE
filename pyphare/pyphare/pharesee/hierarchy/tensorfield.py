from .hierarchy import HierarchyAccessor, PatchHierarchy


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
            hier.ephemerals,
        )

    def __getitem__(self, input):
        if type(input) is HierarchyAccessor:
            return super().__getitem__(input)
        if input in self.__dict__:
            return self.__dict__[input]
        raise IndexError("AnyTensorField.__getitem__ cannot handle input", input)
