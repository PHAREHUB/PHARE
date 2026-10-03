from .hierarchy import PatchHierarchy


class AnyTensorField(PatchHierarchy):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @staticmethod
    def FROM(klass, hier):
        return klass(
            [hier.levels(t) for t in hier.times()],
            hier.domain_box,
            hier.refinement_ratio,
            hier.times(),
            hier.data_files,
            hier.selection_box,
        )

    def __getitem__(self, input):
        if input in self.__dict__:
            return self.__dict__[input]
        raise IndexError("AnyTensorField.__getitem__ cannot handle input", input)
