#
#
#


from pyphare.pharesee.hierarchy import patchlevel


class PatchLevel(patchlevel.PatchLevel):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
