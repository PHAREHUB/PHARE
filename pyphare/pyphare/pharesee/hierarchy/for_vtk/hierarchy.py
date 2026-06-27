#
#
#


from pyphare.pharesee import hierarchy


class PatchHierarchy(hierarchy.PatchHierarchy):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
